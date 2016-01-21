#include <cxxtest/TestSuite.h>
#include "CellBasedSimulationArchiver.hpp"
#include "SmartPointers.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"

#include "DefaultCellProliferativeType.hpp"

#include "CellIdWriter.hpp"
#include "CellAgesWriter.hpp"
#include "VoronoiDataWriter.hpp"
#include "CellMutationStatesWriter.hpp"

#include "ParabolicGrowingDomainPdeModifier.hpp"
#include "MorphogenCellwiseSourceParabolicPde.hpp"
#include "VolumeTrackingModifier.hpp"

#include "MorphogenDependentCellCycleModel.hpp"
#include "CellDataItemWriter.hpp"
#include "OffLatticeSimulation.hpp"
#include "OnLatticeSimulation.hpp"
#include "CellsGenerator.hpp"
#include "RandomCellKiller.hpp"

#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "GeneralisedLinearSpringForce.hpp"

#include "NodeBasedCellPopulation.hpp"
#include "RepulsionForce.hpp"

#include "VertexBasedCellPopulation.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "NagaiHondaForce.hpp"
#include "SimpleTargetAreaModifier.hpp"

#include "PottsBasedCellPopulation.hpp"
#include "PottsMeshGenerator.hpp"
#include "VolumeConstraintPottsUpdateRule.hpp"
#include "AdhesionPottsUpdateRule.hpp"
#include "SurfaceAreaConstraintPottsUpdateRule.hpp"

#include "CaBasedCellPopulation.hpp"
#include "DiffusionCaUpdateRule.hpp"
#include "ShovingCaBasedDivisionRule.hpp"
#include "AdhesionCaSwitchingUpdateRule.hpp"

#include "CommandLineArguments.hpp"
#include "CommandLineArgumentsMocker.hpp"

#include "PetscSetupAndFinalize.hpp"


static const bool M_USING_COMMAND_LINE_ARGS = true;
static const double M_TIME_FOR_SIMULATION = 60; //100
static const double M_NUM_CELLS_ACROSS = 10; // 10
static const double M_UPTAKE_RATE = 0.01; // S in paper
static const double M_DIFFUSION_CONSTANT = 1e-4; // D in paper
static const double M_DUDT_COEFFICIENT = 1.0; // Not used in paper so 1

class TestParabolicMorphogenMonolayers : public AbstractCellBasedWithTimingsTestSuite
{
private:
    void GenerateCells(unsigned num_cells, std::vector<CellPtr>& rCells)
    {
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DefaultCellProliferativeType, p_transit_type);

        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();

        for (unsigned i=0; i<num_cells; i++)
        {
            MorphogenDependentCellCycleModel* p_cycle_model = new MorphogenDependentCellCycleModel();
            p_cycle_model->SetDimension(2);
            p_cycle_model->SetCurrentMass(0.5*(p_gen->ranf()+1.0));
            p_cycle_model->SetMorphogenInfluence(10.0);

            CellPtr p_cell(new Cell(p_state, p_cycle_model));
            p_cell->SetCellProliferativeType(p_transit_type);

            // Note the first few recorded ages will be too short as cells start with some mass.
            p_cell->SetBirthTime(0.0);

            p_cell->InitialiseCellCycleModel();

            // Initial Condition for Morphogen PDE
            p_cell->GetCellData()->SetItem("morphogen",0.0);

            // Set Target Area so dont need to use a growth model in vertex simulations
            p_cell->GetCellData()->SetItem("target area", 1.0);
            rCells.push_back(p_cell);
        }
     }


    /*
     * Simulates diffusion on a growing monolayer.
     */

public:
    void NoTestVertexBasedMonolayer() throw (Exception)
    {
        double sim_index = 0;
        if (M_USING_COMMAND_LINE_ARGS)
        {
            sim_index = (double) atof(CommandLineArguments::Instance()->GetStringCorrespondingToOption("-sim_index").c_str());
        }
        RandomNumberGenerator::Instance()->Reseed(100.0*sim_index);

        //Create output directory
        std::stringstream out;
        out << sim_index;
        std::string output_directory = "ParabolicMonolayers/Vertex/" +  out.str();

        // Create Mesh
        HoneycombVertexMeshGenerator generator(M_NUM_CELLS_ACROSS, M_NUM_CELLS_ACROSS+2);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();
        p_mesh->SetCellRearrangementThreshold(0.1);

        p_mesh->Translate(-(double)M_NUM_CELLS_ACROSS*0.5,-(double)M_NUM_CELLS_ACROSS*0.5);

        // Create Cells
        std::vector<CellPtr> cells;
        GenerateCells(p_mesh->GetNumElements(),cells);

        // Create Population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellAgesWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();
        //Make cell date writer so can pass in variable name
        boost::shared_ptr<CellDataItemWriter<2,2> > p_cell_data_item_writer(new CellDataItemWriter<2,2>("morphogen"));
        cell_population.AddCellWriter(p_cell_data_item_writer);


        // Create Simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory(output_directory);
        simulator.SetDt(1.0/200.0);
        simulator.SetSamplingTimestepMultiple(200);
        simulator.SetEndTime(M_TIME_FOR_SIMULATION);

        simulator.SetOutputDivisionLocations(true);

        // Create Forces and pass to simulation NOTE: these are not the default ones and chosen to give a stable growing monolayer
        MAKE_PTR(NagaiHondaForce<2>, p_force);
        p_force->SetNagaiHondaDeformationEnergyParameter(50.0);
        p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(1.0);
        p_force->SetNagaiHondaCellCellAdhesionEnergyParameter(1.0);
        p_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(10.0);
        simulator.AddForce(p_force);

        // Create Modifiers and pass to simulation

        // Create a pde modifier and pass it to the simulation Add this first so in place for SimpleTargetArea one (calls cell pop update)

		// Make the Pde and BCS
        MorphogenCellwiseSourceParabolicPde<2> pde(cell_population, M_DUDT_COEFFICIENT,M_DIFFUSION_CONSTANT,M_UPTAKE_RATE);
		ConstBoundaryCondition<2> bc(0.0);
		ParabolicPdeAndBoundaryConditions<2> pde_and_bc(&pde, &bc, true);
		pde_and_bc.SetDependentVariableName("morphogen");

		// Create a PDE Modifier object using this pde and bcs object
		MAKE_PTR_ARGS(ParabolicGrowingDomainPdeModifier<2>, p_pde_modifier, (&pde_and_bc));
		simulator.AddSimulationModifier(p_pde_modifier);

        simulator.Solve();
    }

    void NoTestNodeBasedMonolayer() throw (Exception)
    {
        double sim_index = 0;
        if (M_USING_COMMAND_LINE_ARGS)
        {
            sim_index = (double) atof(CommandLineArguments::Instance()->GetStringCorrespondingToOption("-sim_index").c_str());
        }
        RandomNumberGenerator::Instance()->Reseed(100.0*sim_index);

        double cut_off_length = 1.5; //this is the default

        //Create output directory
        std::stringstream out;
        out << sim_index;
        std::string output_directory = "ParabolicMonolayers/Node/" +  out.str();

        HoneycombMeshGenerator generator(M_NUM_CELLS_ACROSS, M_NUM_CELLS_ACROSS+2,0);
        MutableMesh<2,2>* p_generating_mesh = generator.GetMesh();
        NodesOnlyMesh<2>* p_mesh = new NodesOnlyMesh<2>;
        p_mesh->ConstructNodesWithoutMesh(*p_generating_mesh, cut_off_length);

        p_mesh->Translate(-(double)M_NUM_CELLS_ACROSS*0.5,-(double)M_NUM_CELLS_ACROSS*0.5);

        std::vector<CellPtr> cells;
        GenerateCells(p_mesh->GetNumNodes(),cells);

        NodeBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellAgesWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();
        //Make cell date writer so can pass in variable name
        boost::shared_ptr<CellDataItemWriter<2,2> > p_cell_data_item_writer(new CellDataItemWriter<2,2>("morphogen"));
        cell_population.AddCellWriter(p_cell_data_item_writer);

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory(output_directory);
        simulator.SetDt(1.0/200.0);
        simulator.SetSamplingTimestepMultiple(200);
        simulator.SetEndTime(M_TIME_FOR_SIMULATION);

        simulator.SetOutputDivisionLocations(true);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetMeinekeSpringStiffness(50.0);
        p_linear_force->SetCutOffLength(cut_off_length);
        simulator.AddForce(p_linear_force);

        // Make the Pde and BCS
        MorphogenCellwiseSourceParabolicPde<2> pde(cell_population, M_DUDT_COEFFICIENT,M_DIFFUSION_CONSTANT,M_UPTAKE_RATE);
		ConstBoundaryCondition<2> bc(0.0);
		ParabolicPdeAndBoundaryConditions<2> pde_and_bc(&pde, &bc, true);
		pde_and_bc.SetDependentVariableName("morphogen");

		// Create a PDE Modifier object using this pde and bcs object
		MAKE_PTR_ARGS(ParabolicGrowingDomainPdeModifier<2>, p_pde_modifier, (&pde_and_bc));
		simulator.AddSimulationModifier(p_pde_modifier);

        simulator.Solve();

        delete p_mesh; // to stop memory leaks
    }

    void NoTestMeshBasedMonolayer() throw (Exception)
    {
        double sim_index = 0;
        if (M_USING_COMMAND_LINE_ARGS)
        {
            sim_index = (double) atof(CommandLineArguments::Instance()->GetStringCorrespondingToOption("-sim_index").c_str());
        }
        RandomNumberGenerator::Instance()->Reseed(100.0*sim_index);

        //Create output directory
        std::stringstream out;
        out << sim_index;
        std::string output_directory = "ParabolicMonolayers/MeshPoint/" +  out.str();

        HoneycombMeshGenerator generator(M_NUM_CELLS_ACROSS,M_NUM_CELLS_ACROSS+2,0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        p_mesh->Translate(-(double)M_NUM_CELLS_ACROSS*0.5,-(double)M_NUM_CELLS_ACROSS*0.5);

        std::vector<CellPtr> cells;
        GenerateCells(p_mesh->GetNumNodes(),cells);

        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set population to output all data to results files
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellAgesWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();
        //Make cell date writer so can pass in variabl name
        boost::shared_ptr<CellDataItemWriter<2,2> > p_cell_data_item_writer(new CellDataItemWriter<2,2>("morphogen"));
        cell_population.AddCellWriter(p_cell_data_item_writer);

        cell_population.SetWriteVtkAsPoints(true);
        cell_population.AddPopulationWriter<VoronoiDataWriter>();


        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory(output_directory);
        simulator.SetDt(1.0/200.0);
        simulator.SetSamplingTimestepMultiple(200);
        simulator.SetEndTime(M_TIME_FOR_SIMULATION);

        simulator.SetOutputDivisionLocations(true);

        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetMeinekeSpringStiffness(50.0);
        p_linear_force->SetCutOffLength(1.5);
        simulator.AddForce(p_linear_force);

        // Make the Pde and BCS
        MorphogenCellwiseSourceParabolicPde<2> pde(cell_population, M_DUDT_COEFFICIENT,M_DIFFUSION_CONSTANT,M_UPTAKE_RATE);
		ConstBoundaryCondition<2> bc(0.0);
		ParabolicPdeAndBoundaryConditions<2> pde_and_bc(&pde, &bc, true);
		pde_and_bc.SetDependentVariableName("morphogen");

		// Create a PDE Modifier object using this pde and bcs object
		MAKE_PTR_ARGS(ParabolicGrowingDomainPdeModifier<2>, p_pde_modifier, (&pde_and_bc));
		simulator.AddSimulationModifier(p_pde_modifier);

        simulator.Solve();
    }

    void NoTestPottsBasedMonolayer() throw (Exception)
    {
        double sim_index = 0;
        if (M_USING_COMMAND_LINE_ARGS)
        {
            sim_index = (double) atof(CommandLineArguments::Instance()->GetStringCorrespondingToOption("-sim_index").c_str());
        }
        RandomNumberGenerator::Instance()->Reseed(100.0*sim_index);

        //Create output directory
        std::stringstream out;
        out << sim_index;
        std::string output_directory = "ParabolicMonolayers/Potts/" +  out.str();

        unsigned cell_width = 4;
        unsigned domain_width = M_NUM_CELLS_ACROSS*cell_width*20;
        PottsMeshGenerator<2> generator(domain_width, M_NUM_CELLS_ACROSS, cell_width, domain_width, M_NUM_CELLS_ACROSS, cell_width);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        p_mesh->Translate(-(double)domain_width*0.5+0.5,-(double)domain_width*0.5+0.5);

        p_mesh->Scale(0.25,0.25);

        std::vector<CellPtr> cells;
        GenerateCells(p_mesh->GetNumElements(),cells);

        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetTemperature(0.1);
        cell_population.SetNumSweepsPerTimestep(1);

        // Set population to output all data to results files
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellAgesWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();
        //Make cell date writer so can pass in variable name
        boost::shared_ptr<CellDataItemWriter<2,2> > p_cell_data_item_writer(new CellDataItemWriter<2,2>("morphogen"));
        cell_population.AddCellWriter(p_cell_data_item_writer);

        // Set the Temperature
        cell_population.SetTemperature(0.1); //Default is 0.1

        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory(output_directory);
        simulator.SetDt(1.0/200.0);
        simulator.SetSamplingTimestepMultiple(200);
        simulator.SetEndTime(M_TIME_FOR_SIMULATION);

        simulator.SetOutputDivisionLocations(true);

        // Create update rules and pass to the simulation
        MAKE_PTR(VolumeConstraintPottsUpdateRule<2>, p_volume_constraint_update_rule);
        p_volume_constraint_update_rule->SetMatureCellTargetVolume(16); // i.e 4x4 cells
        p_volume_constraint_update_rule->SetDeformationEnergyParameter(0.1);
        simulator.AddPottsUpdateRule(p_volume_constraint_update_rule);

        MAKE_PTR(SurfaceAreaConstraintPottsUpdateRule<2>, p_surface_constraint_update_rule);
        p_surface_constraint_update_rule->SetMatureCellTargetSurfaceArea(16); // i.e 4x4 cells
        p_surface_constraint_update_rule->SetDeformationEnergyParameter(0.01);
        simulator.AddPottsUpdateRule(p_surface_constraint_update_rule);

        MAKE_PTR(AdhesionPottsUpdateRule<2>, p_adhesion_update_rule);
        p_adhesion_update_rule->SetCellCellAdhesionEnergyParameter(0.1);
        p_adhesion_update_rule->SetCellBoundaryAdhesionEnergyParameter(0.2);
        simulator.AddPottsUpdateRule(p_adhesion_update_rule);

        // Make the Pde and BCS
        //MorphogenCellwiseSourceParabolicPde<2> pde(cell_population, M_DUDT_COEFFICIENT,(double)cell_width*(double)cell_width*M_DIFFUSION_CONSTANT,M_UPTAKE_RATE, 8.0);
        MorphogenCellwiseSourceParabolicPde<2> pde(cell_population, M_DUDT_COEFFICIENT,M_DIFFUSION_CONSTANT,M_UPTAKE_RATE);
        ConstBoundaryCondition<2> bc(0.0);
		ParabolicPdeAndBoundaryConditions<2> pde_and_bc(&pde, &bc, true);
		pde_and_bc.SetDependentVariableName("morphogen");

		// Create a PDE Modifier object using this pde and bcs object
		MAKE_PTR_ARGS(ParabolicGrowingDomainPdeModifier<2>, p_pde_modifier, (&pde_and_bc));
		simulator.AddSimulationModifier(p_pde_modifier);

        simulator.Solve();
    }

    void TestCaBasedMonolayer() throw (Exception)
    {
        double sim_index = 0;
        if (M_USING_COMMAND_LINE_ARGS)
        {
            sim_index = (double) atof(CommandLineArguments::Instance()->GetStringCorrespondingToOption("-sim_index").c_str());
        }
        RandomNumberGenerator::Instance()->Reseed(100.0*sim_index);

        //Create output directory
        std::stringstream out;
        out << sim_index;
        std::string output_directory = "ParabolicMonolayers/Ca/" +  out.str();

        // Create a simple 2D PottsMesh
        unsigned domain_wide = 10*M_NUM_CELLS_ACROSS;

        PottsMeshGenerator<2> generator(domain_wide, 0, 0, domain_wide, 0, 0);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        p_mesh->Translate(-(double)domain_wide*0.5 + 0.5,-(double)domain_wide*0.5 + 0.5);

        // Specify where cells lie
        std::vector<unsigned> location_indices;
        for (unsigned i=0; i<M_NUM_CELLS_ACROSS; i++)
        {
            for (unsigned j=0; j<M_NUM_CELLS_ACROSS; j++)
            {
                unsigned offset = (domain_wide+1) * (domain_wide-M_NUM_CELLS_ACROSS)/2;
                location_indices.push_back(offset + j + i * domain_wide );
            }
        }

        std::vector<CellPtr> cells;
        GenerateCells(location_indices.size(),cells);

        // Create cell population
        CaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);

        // Set population to output all data to results files
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellAgesWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();
        //Make cell date writer so can pass in variable name
        boost::shared_ptr<CellDataItemWriter<2,2> > p_cell_data_item_writer(new CellDataItemWriter<2,2>("morphogen"));
        cell_population.AddCellWriter(p_cell_data_item_writer);

        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory(output_directory);
        simulator.SetDt(1/200.0);
        simulator.SetSamplingTimestepMultiple(200);
        simulator.SetEndTime(M_TIME_FOR_SIMULATION);

        simulator.SetOutputDivisionLocations(true);

        // Add Division Rule
        boost::shared_ptr<AbstractCaBasedDivisionRule<2> > p_division_rule(new ShovingCaBasedDivisionRule<2>());
        cell_population.SetCaBasedDivisionRule(p_division_rule);

        // Add switching Update Rule to smooth out the edge of the monolayer note no Temp as don't want random switches
        MAKE_PTR(AdhesionCaSwitchingUpdateRule<2u>, p_switching_update_rule);
        p_switching_update_rule->SetCellCellAdhesionEnergyParameter(0.1);
        p_switching_update_rule->SetCellBoundaryAdhesionEnergyParameter(0.2);
        p_switching_update_rule->SetTemperature(0.0); //
        simulator.AddCaSwitchingUpdateRule(p_switching_update_rule);

        // Make the Pde and BCS
        MorphogenCellwiseSourceParabolicPde<2> pde(cell_population, M_DUDT_COEFFICIENT,M_DIFFUSION_CONSTANT,M_UPTAKE_RATE);
        ConstBoundaryCondition<2> bc(0.0);
        ParabolicPdeAndBoundaryConditions<2> pde_and_bc(&pde, &bc, true);
        pde_and_bc.SetDependentVariableName("morphogen");

        // Create a PDE Modifier object using this pde and bcs object
        MAKE_PTR_ARGS(ParabolicGrowingDomainPdeModifier<2>, p_pde_modifier, (&pde_and_bc));
        simulator.AddSimulationModifier(p_pde_modifier);

        simulator.Solve();
    }
};
