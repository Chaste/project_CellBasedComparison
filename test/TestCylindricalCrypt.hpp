#ifndef TESTCRYPTSIMULATIONCOMPARISON_HPP_
#define TESTCRYPTSIMULATIONCOMPARISON_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before any other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "CylindricalHoneycombVertexMeshGenerator.hpp"
#include "CylindricalHoneycombMeshGenerator.hpp"
#include "PottsMeshGenerator.hpp"
#include "Cylindrical2dNodesOnlyMesh.hpp"

#include "CellsGenerator.hpp"

#include "SimpleWntContactInhibitionCellCycleModel.hpp"
#include "WntConcentration.hpp"

#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"

#include "CellProliferativeTypesCountWriter.hpp"
#include "CellIdWriter.hpp"
#include "CellVolumesWriter.hpp"
#include "CellAncestorWriter.hpp"

#include "OffLatticeSimulation.hpp"
#include "OnLatticeSimulation.hpp"

#include "NagaiHondaForce.hpp"
#include "RepulsionForce.hpp"
#include "DiffusionCaUpdateRule.hpp"
#include "VolumeConstraintPottsUpdateRule.hpp"
#include "AdhesionPottsUpdateRule.hpp"
#include "SurfaceAreaConstraintPottsUpdateRule.hpp"


#include "SimpleTargetAreaModifier.hpp"
#include "VolumeTrackingModifier.hpp"

#include "PlaneBasedCellKiller.hpp"

#include "PlaneBoundaryCondition.hpp"

#include "CryptShovingCaBasedDivisionRule.hpp"

#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "Warnings.hpp"

static const bool M_USING_COMMAND_LINE_ARGS = true;
static const double M_END_STEADY_STATE = 100.0; //100
static const double M_END_TIME = 1100; //1100
static const double M_CRYPT_DIAMETER = 16;
static const double M_CRYPT_LENGTH = 12;

class TestCryptSimulationComparison : public AbstractCellBasedWithTimingsTestSuite
{
private:

    void GenerateCells(unsigned num_cells, std::vector<CellPtr>& rCells, double equilibriumVolume, double quiescentVolumeFraction)
    {
        double typical_cell_cycle_duration = 12.0;

        boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_cell_type(CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>());

        for (unsigned i=0; i<num_cells; i++)
        {
            SimpleWntContactInhibitionCellCycleModel* p_model = new SimpleWntContactInhibitionCellCycleModel();
            p_model->SetDimension(2);
            p_model->SetEquilibriumVolume(equilibriumVolume);
            p_model->SetQuiescentVolumeFraction(quiescentVolumeFraction);
            p_model->SetWntThreshold(0.5); //0.8 -> CI // 0.1 -> No CI!!!!


            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetCellProliferativeType(p_cell_type);
            double birth_time = - RandomNumberGenerator::Instance()->ranf() * typical_cell_cycle_duration;
            p_cell->SetBirthTime(birth_time);

            // Set Target Area so dont need to use a growth model in vertex simulations
            p_cell->GetCellData()->SetItem("target area", 1.0);

            rCells.push_back(p_cell);
        }
    }

public:

    void NoTestVertexCrypt() throw (Exception)
    {
        double sim_index = 0;
        double contact_inhibition_level = 0.8;
        if (M_USING_COMMAND_LINE_ARGS)
        {
          sim_index = (double) atof(CommandLineArguments::Instance()->GetStringCorrespondingToOption("-sim_index").c_str());
          contact_inhibition_level = (double) atof(CommandLineArguments::Instance()->GetStringCorrespondingToOption("-CI").c_str());
        }
        RandomNumberGenerator::Instance()->Reseed(100.0*sim_index);

        //Create output directory
        std::stringstream out;
        out << sim_index << "_CI_" << contact_inhibition_level;
        std::string output_directory = "CylindricalCrypt/Vertex/" +  out.str();

        // Create mesh
        CylindricalHoneycombVertexMeshGenerator generator(M_CRYPT_DIAMETER, M_CRYPT_LENGTH, true);
        Cylindrical2dVertexMesh* p_mesh = generator.GetCylindricalMesh();
        p_mesh->SetCellRearrangementThreshold(0.1);


        // Create cells
        std::vector<CellPtr> cells;
        GenerateCells(p_mesh->GetNumElements(),cells,1.0,contact_inhibition_level); //mature_volume = 1.0

        // Create tissue
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
        cell_population.AddCellWriter<CellVolumesWriter>();
        cell_population.AddCellWriter<CellIdWriter>();

        // Create an instance of a Wnt concentration
        WntConcentration<2>::Instance()->SetType(LINEAR);
        WntConcentration<2>::Instance()->SetCellPopulation(cell_population);
        WntConcentration<2>::Instance()->SetCryptLength(M_CRYPT_LENGTH);

        // Create crypt simulation from cell population
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetDt(1.0/200.0);
        simulator.SetSamplingTimestepMultiple(200);
        simulator.SetEndTime(M_END_STEADY_STATE);
        simulator.SetOutputDirectory(output_directory);
        simulator.SetOutputDivisionLocations(true);
        simulator.SetOutputCellVelocities(true);
        cell_population.AddCellWriter<CellAncestorWriter>();

        // Add volume tracking modifier
        MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        // Create Forces and pass to simulation NOTE : these are not the default ones and chosen to give a stable growing monolayer
        MAKE_PTR(NagaiHondaForce<2>, p_force);
        p_force->SetNagaiHondaDeformationEnergyParameter(50.0);
        p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(1.0);
        p_force->SetNagaiHondaCellCellAdhesionEnergyParameter(1.0);
        p_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(10.0);
        simulator.AddForce(p_force);

        // A NagaiHondaForce has to be used together with an AbstractTargetAreaModifier
//        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
//        simulator.AddSimulationModifier(p_growth_modifier);

        // Solid base Boundary condition
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bcs, (&cell_population, zero_vector<double>(2), -unit_vector<double>(2,1)));
        p_bcs->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(p_bcs);

        // Sloughing killer
        MAKE_PTR_ARGS(PlaneBasedCellKiller<2>, p_killer, (&cell_population, M_CRYPT_LENGTH*unit_vector<double>(2,1), unit_vector<double>(2,1)));
        simulator.AddCellKiller(p_killer);

        // Run simulation
        simulator.Solve();

        // Mark Ancestors
        simulator.SetEndTime(M_END_TIME);
        simulator.rGetCellPopulation().SetCellAncestorsToLocationIndices();

        // Run simulation to new end time
        simulator.Solve();

        // Clear singletons
        WntConcentration<2>::Instance()->Destroy();
        Warnings::Instance()->QuietDestroy();
    }

    void NoTestMeshBasedCrypt() throw (Exception)
    {
        double sim_index = 0;
        double contact_inhibition_level = 0.8;
        if (M_USING_COMMAND_LINE_ARGS)
        {
          sim_index = (double) atof(CommandLineArguments::Instance()->GetStringCorrespondingToOption("-sim_index").c_str());
          contact_inhibition_level = (double) atof(CommandLineArguments::Instance()->GetStringCorrespondingToOption("-CI").c_str());
        }
        RandomNumberGenerator::Instance()->Reseed(100.0*sim_index);
        //Create output directory
        std::stringstream out;
        out << sim_index << "_CI_" << contact_inhibition_level;
        std::string output_directory = "CylindricalCrypt/Mesh/" +  out.str();

        // Create mesh
        unsigned thickness_of_ghost_layer = 2;

        CylindricalHoneycombMeshGenerator generator(M_CRYPT_DIAMETER, M_CRYPT_LENGTH, thickness_of_ghost_layer);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        // Get location indices corresponding to real cells
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Create cells
        std::vector<CellPtr> cells;
        GenerateCells(location_indices.size(),cells,sqrt(3.0)/2.0,contact_inhibition_level);  //mature_volume = sqrt(3.0)/2.0

        // Create tissue
        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);
        cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
        cell_population.AddCellWriter<CellVolumesWriter>();
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellAncestorWriter>();

        // Create an instance of a Wnt concentration
        WntConcentration<2>::Instance()->SetType(LINEAR);
        WntConcentration<2>::Instance()->SetCellPopulation(cell_population);
        WntConcentration<2>::Instance()->SetCryptLength(M_CRYPT_LENGTH);

        // Create simulation from cell population
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetDt(1.0/200.0);
        simulator.SetSamplingTimestepMultiple(200);
        simulator.SetEndTime(M_END_STEADY_STATE);
        simulator.SetOutputDirectory(output_directory);
        simulator.SetOutputDivisionLocations(true);
        simulator.SetOutputCellVelocities(true);

        // Add volume tracking Modifier
        MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetMeinekeSpringStiffness(50.0);
        simulator.AddForce(p_linear_force);

        // Solid base boundary condition
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bcs, (&cell_population, zero_vector<double>(2), -unit_vector<double>(2,1)));
        p_bcs->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(p_bcs);

        // Sloughing killer
        MAKE_PTR_ARGS(PlaneBasedCellKiller<2>, p_killer, (&cell_population, (M_CRYPT_LENGTH-0.5)*unit_vector<double>(2,1), unit_vector<double>(2,1)));
        simulator.AddCellKiller(p_killer);

        // Run simulation
        simulator.Solve();

        // Mark Ancestors
        simulator.SetEndTime(M_END_TIME);
        simulator.rGetCellPopulation().SetCellAncestorsToLocationIndices();

        // Run simulation to new end time
        simulator.Solve();

        // Clear singletons
        WntConcentration<2>::Instance()->Destroy();
    }

    void NoTestNodeBasedCrypt() throw (Exception)
    {
        double sim_index = 0;
        double contact_inhibition_level = 0.8;
        if (M_USING_COMMAND_LINE_ARGS)
        {
          sim_index = (double) atof(CommandLineArguments::Instance()->GetStringCorrespondingToOption("-sim_index").c_str());
          contact_inhibition_level = (double) atof(CommandLineArguments::Instance()->GetStringCorrespondingToOption("-CI").c_str());
        }
        RandomNumberGenerator::Instance()->Reseed(100.0*sim_index);

        double cut_off_length = 1.5; //this is the default

        //Create output directory
        std::stringstream out;
        out << sim_index << "_CI_" << contact_inhibition_level;
        std::string output_directory = "CylindricalCrypt/Node/" +  out.str();

        // Create a simple mesh
        HoneycombMeshGenerator generator(M_CRYPT_DIAMETER, M_CRYPT_LENGTH, 0);
        TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();

        // Convert this to a Cylindrical2dNodesOnlyMesh
        Cylindrical2dNodesOnlyMesh* p_mesh = new Cylindrical2dNodesOnlyMesh(M_CRYPT_DIAMETER);
        p_mesh->ConstructNodesWithoutMesh(*p_generating_mesh,M_CRYPT_DIAMETER);

        // Create cells
        std::vector<CellPtr> cells;
        GenerateCells(p_mesh->GetNumNodes(), cells, M_PI*0.25,contact_inhibition_level); // mature volume: M_PI*0.25 as r=0.5

        // Create a node-based cell population
        NodeBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
        cell_population.AddCellWriter<CellVolumesWriter>();
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellAncestorWriter>();

        for (unsigned index = 0; index < cell_population.rGetMesh().GetNumNodes(); index++)
        {
            cell_population.rGetMesh().GetNode(index)->SetRadius(0.5);
        }

        // Create an instance of a Wnt concentration
        WntConcentration<2>::Instance()->SetType(LINEAR);
        WntConcentration<2>::Instance()->SetCellPopulation(cell_population);
        WntConcentration<2>::Instance()->SetCryptLength(M_CRYPT_LENGTH);

        // Create simulation from cell population
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetDt(1.0/200.0);
        simulator.SetSamplingTimestepMultiple(200);
        simulator.SetEndTime(M_END_STEADY_STATE);
        simulator.SetOutputDirectory(output_directory);
        simulator.SetOutputDivisionLocations(true);
        simulator.SetOutputCellVelocities(true);

        // Add volume tracking modifier
        MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetMeinekeSpringStiffness(50.0);
        p_linear_force->SetCutOffLength(cut_off_length);
        simulator.AddForce(p_linear_force);

        // Solid base boundary condition
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bcs, (&cell_population, zero_vector<double>(2), -unit_vector<double>(2,1)));
        p_bcs->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(p_bcs);

        // Sloughing killer
        MAKE_PTR_ARGS(PlaneBasedCellKiller<2>, p_killer, (&cell_population, (M_CRYPT_LENGTH-0.5)*unit_vector<double>(2,1), unit_vector<double>(2,1)));
        simulator.AddCellKiller(p_killer);

        // Run simulation
        simulator.Solve();

        // Mark Ancestors
        simulator.SetEndTime(M_END_TIME);
        simulator.rGetCellPopulation().SetCellAncestorsToLocationIndices();

        // Run simulation to new end time
        simulator.Solve();

        // Clear memory
        delete p_mesh;

        // Clear singletons
        WntConcentration<2>::Instance()->Destroy();
    }

    void TestPottsCrypt() throw (Exception)
    {
        double sim_index = 0;
        double contact_inhibition_level = 0.8;
        if (M_USING_COMMAND_LINE_ARGS)
        {
          sim_index = (double) atof(CommandLineArguments::Instance()->GetStringCorrespondingToOption("-sim_index").c_str());
          contact_inhibition_level = (double) atof(CommandLineArguments::Instance()->GetStringCorrespondingToOption("-CI").c_str());
        }
        RandomNumberGenerator::Instance()->Reseed(100.0*sim_index);

        //Create output directory
        std::stringstream out;
        out << sim_index << "_CI_" << contact_inhibition_level;
        std::string output_directory = "CylindricalCrypt/Potts/" +  out.str();

        unsigned cell_width = 4;

        // Create a simple 2D PottsMesh (periodic in x)
        PottsMeshGenerator<2> generator( M_CRYPT_DIAMETER*cell_width, M_CRYPT_DIAMETER, cell_width, (M_CRYPT_LENGTH+2)*cell_width, M_CRYPT_LENGTH, cell_width, 1, 1, 1, true, true); //Dtart from bottom left and periodic
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        GenerateCells(p_mesh->GetNumElements(),cells,cell_width*cell_width,contact_inhibition_level); // mature volume = 16.0 LSs

        // Create cell population
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetTemperature(0.1);
        cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
        cell_population.AddCellWriter<CellVolumesWriter>();
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellAncestorWriter>();

        // Set the Temperature
        cell_population.SetTemperature(0.1); //Default is 0.1


        // Create an instance of a Wnt concentration
        WntConcentration<2>::Instance()->SetType(LINEAR);
        WntConcentration<2>::Instance()->SetCellPopulation(cell_population);
        WntConcentration<2>::Instance()->SetCryptLength(M_CRYPT_LENGTH*cell_width);

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory(output_directory);
        simulator.SetDt(0.01);
        simulator.SetSamplingTimestepMultiple(100);
        simulator.SetEndTime(M_END_STEADY_STATE);
        simulator.SetOutputDivisionLocations(true);
        simulator.SetOutputCellVelocities(true);

        // Add volume tracking modifier
        MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        // Sloughing killer
        MAKE_PTR_ARGS(PlaneBasedCellKiller<2>, p_killer, (&cell_population, cell_width*M_CRYPT_LENGTH*unit_vector<double>(2,1), unit_vector<double>(2,1)));
        simulator.AddCellKiller(p_killer);

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

        // Run simulation
        simulator.Solve();

        // Mark Ancestors
        simulator.SetEndTime(M_END_TIME);
        simulator.rGetCellPopulation().SetCellAncestorsToLocationIndices();

        // Run simulation to new end time
        simulator.Solve();

        // Clear singletons
        WntConcentration<2>::Instance()->Destroy();
    }

    void NoTestCaCrypt() throw (Exception)
    {
        double sim_index = 0;
        double contact_inhibition_level = 0.8;
        if (M_USING_COMMAND_LINE_ARGS)
        {
          sim_index = (double) atof(CommandLineArguments::Instance()->GetStringCorrespondingToOption("-sim_index").c_str());
          contact_inhibition_level = (double) atof(CommandLineArguments::Instance()->GetStringCorrespondingToOption("-CI").c_str());
        }
        RandomNumberGenerator::Instance()->Reseed(100.0*sim_index);

        //Create output directory
        std::stringstream out;
        out << sim_index << "_CI_" << contact_inhibition_level;
        std::string output_directory = "CylindricalCrypt/Ca/" +  out.str();

        // Create a simple 2D PottsMesh (periodic in x)
        PottsMeshGenerator<2> generator(M_CRYPT_DIAMETER, 0, 0, M_CRYPT_LENGTH*3, 0, 0, 1, 0, 0, false, true);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Specify where cells lie
        std::vector<unsigned> location_indices;
        for (unsigned index=0; index<(unsigned)M_CRYPT_DIAMETER*(unsigned)M_CRYPT_LENGTH; index++)
        {
            location_indices.push_back(index);
        }

//        // Create cells
        std::vector<CellPtr> cells;
        GenerateCells(location_indices.size(),cells,1.0,contact_inhibition_level); //Mature volume = 1 LS

        // Create cell population
        CaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);
        cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
        cell_population.AddCellWriter<CellVolumesWriter>();
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellAncestorWriter>();

        // Create an instance of a Wnt concentration
        WntConcentration<2>::Instance()->SetType(LINEAR);
        WntConcentration<2>::Instance()->SetCellPopulation(cell_population);
        WntConcentration<2>::Instance()->SetCryptLength(M_CRYPT_LENGTH);

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory(output_directory);
        simulator.SetDt(0.01);
        simulator.SetSamplingTimestepMultiple(100);
        simulator.SetEndTime(M_END_STEADY_STATE);
        simulator.SetOutputDivisionLocations(true);
        simulator.SetOutputCellVelocities(true);

        // Add Volume tracking modifier
        MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        // Add Division Rule
        boost::shared_ptr<AbstractCaBasedDivisionRule<2> > p_division_rule(new CryptShovingCaBasedDivisionRule());
        cell_population.SetCaBasedDivisionRule(p_division_rule);


        // Sloughing killer
        MAKE_PTR_ARGS(PlaneBasedCellKiller<2>, p_killer, (&cell_population, M_CRYPT_LENGTH*unit_vector<double>(2,1), unit_vector<double>(2,1)));
        simulator.AddCellKiller(p_killer);

        // Run simulation
        simulator.Solve();

        // Mark Ancestors
        simulator.SetEndTime(M_END_TIME);
        simulator.rGetCellPopulation().SetCellAncestorsToLocationIndices();

        // Run simulation to new end time
        simulator.Solve();

        // Clear singletons
        WntConcentration<2>::Instance()->Destroy();
    }
};

#endif /*TESTCRYPTSIMULATIONCOMPARISON_HPP_*/
