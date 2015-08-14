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
#include "CellMorphogenWriter.hpp"
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

#include "CommandLineArguments.hpp"
#include "CommandLineArgumentsMocker.hpp"

#include "PetscSetupAndFinalize.hpp"


static const bool M_USING_COMMAND_LINE_ARGS = true;
static const double M_TIME_FOR_SIMULATION = 1.0;
static const double M_NUM_CELLS_ACROSS = 10; // this ^2 cells
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

            // TODO REMOVE THIS AS MAY SKEW AGE PLOTS!
            double typical_cell_cycle_duration = 12.0;
            double birth_time = - p_gen->ranf() * typical_cell_cycle_duration;
            p_cell->SetBirthTime(birth_time);

            p_cell->InitialiseCellCycleModel();
            p_cell->GetCellData()->SetItem("morphogen",0.0);

            rCells.push_back(p_cell);
        }
     }


    /*
     * Simulates diffusion on a growing monolayer.
     */

public:
    void TestVertexBasedMonolayer() throw (Exception)
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

        p_mesh->Translate(-(double)M_NUM_CELLS_ACROSS*0.5,-(double)M_NUM_CELLS_ACROSS*0.5);

        // Create Cells
        std::vector<CellPtr> cells;
        GenerateCells(p_mesh->GetNumElements(),cells);

        // Create Population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellAgesWriter>();
        cell_population.AddCellWriter<CellMorphogenWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();

        // Create Simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory(output_directory);
        simulator.SetDt(1.0/200.0);
        simulator.SetSamplingTimestepMultiple(200);
        simulator.SetEndTime(M_TIME_FOR_SIMULATION);//20

        simulator.SetOutputDivisionLocations(true);

        // Create Forces and pass to simulation NOTE: PARAMETERS CHOSEN TO GET CIRCULAR MONOLAYER
        MAKE_PTR(NagaiHondaForce<2>, p_force);
        p_force->SetNagaiHondaDeformationEnergyParameter(55.0);
        p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(0.0);
        p_force->SetNagaiHondaCellCellAdhesionEnergyParameter(6.0);
        p_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(12.0);
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

        // A NagaiHondaForce has to be used together with an AbstractTargetAreaModifier
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        simulator.Solve();
    }

    void TestNodeBasedMonolayer() throw (Exception)
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
        std::string output_directory = "ParabolicMonolayers/Node/" +  out.str();

        HoneycombMeshGenerator generator(M_NUM_CELLS_ACROSS, M_NUM_CELLS_ACROSS+2,0);
        MutableMesh<2,2>* p_generating_mesh = generator.GetMesh();
        NodesOnlyMesh<2>* p_mesh = new NodesOnlyMesh<2>;
        p_mesh->ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);

        p_mesh->Translate(-(double)M_NUM_CELLS_ACROSS*0.5,-(double)M_NUM_CELLS_ACROSS*0.5);

        std::vector<CellPtr> cells;
        GenerateCells(p_mesh->GetNumNodes(),cells);

        NodeBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellAgesWriter>();
        cell_population.AddCellWriter<CellMorphogenWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory(output_directory);
        simulator.SetDt(1.0/120.0);
        simulator.SetSamplingTimestepMultiple(120);
        simulator.SetEndTime(M_TIME_FOR_SIMULATION);

        simulator.SetOutputDivisionLocations(true);

        MAKE_PTR(RepulsionForce<2>, p_force);
        simulator.AddForce(p_force);

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

    void TestMeshBasedMonolayer() throw (Exception)
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
        cell_population.AddCellWriter<CellMorphogenWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();

        cell_population.SetWriteVtkAsPoints(true);
        cell_population.AddPopulationWriter<VoronoiDataWriter>();


        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory(output_directory);
        simulator.SetDt(1.0/120.0);
        simulator.SetSamplingTimestepMultiple(120);
        simulator.SetEndTime(M_TIME_FOR_SIMULATION);

        simulator.SetOutputDivisionLocations(true);

        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        p_force->SetCutOffLength(1.5);
        simulator.AddForce(p_force);

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

    void TestPottsBasedMonolayer() throw (Exception)
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
        unsigned domain_width = M_NUM_CELLS_ACROSS*cell_width*5;
        PottsMeshGenerator<2> generator(domain_width, M_NUM_CELLS_ACROSS, cell_width, domain_width, M_NUM_CELLS_ACROSS, cell_width);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        p_mesh->Translate(-(double)domain_width*0.5+0.5,-(double)domain_width*0.5+0.5);

        std::vector<CellPtr> cells;
        GenerateCells(p_mesh->GetNumElements(),cells);

        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetTemperature(0.1);
        cell_population.SetNumSweepsPerTimestep(100);

        // Set population to output all data to results files
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellAgesWriter>();
        cell_population.AddCellWriter<CellMorphogenWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();

        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory(output_directory);
        simulator.SetDt(1.0);
        simulator.SetSamplingTimestepMultiple(1);
        simulator.SetEndTime(M_TIME_FOR_SIMULATION);

        simulator.SetOutputDivisionLocations(true);

        MAKE_PTR(VolumeConstraintPottsUpdateRule<2>, p_volume_constraint_update_rule);
        simulator.AddPottsUpdateRule(p_volume_constraint_update_rule);
        MAKE_PTR(SurfaceAreaConstraintPottsUpdateRule<2>, p_surface_area_update_rule);
        simulator.AddPottsUpdateRule(p_surface_area_update_rule);
        MAKE_PTR(AdhesionPottsUpdateRule<2>, p_adhesion_update_rule);
        simulator.AddPottsUpdateRule(p_adhesion_update_rule);

        // Make the Pde and BCS
        MorphogenCellwiseSourceParabolicPde<2> pde(cell_population, M_DUDT_COEFFICIENT,(double)cell_width*(double)cell_width*M_DIFFUSION_CONSTANT,M_UPTAKE_RATE, 8.0);
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
        unsigned domain_wide = 5*M_NUM_CELLS_ACROSS;

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
        cell_population.AddCellWriter<CellMorphogenWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();

        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory(output_directory);
        simulator.SetDt(1.0);
        simulator.SetSamplingTimestepMultiple(1);
        simulator.SetEndTime(M_TIME_FOR_SIMULATION);

        simulator.SetOutputDivisionLocations(true);

        // Add Division Rule
        boost::shared_ptr<AbstractCaBasedDivisionRule<2> > p_division_rule(new ShovingCaBasedDivisionRule<2>());
        cell_population.SetCaBasedDivisionRule(p_division_rule);

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
