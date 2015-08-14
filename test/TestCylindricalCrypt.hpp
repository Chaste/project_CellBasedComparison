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

#include "ContactInhibitionCellCycleModel.hpp"
#include "StochasticDurationCellCycleModel.hpp"

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

#include "SimpleTargetAreaModifier.hpp"
#include "VolumeTrackingModifier.hpp"

#include "PlaneBasedCellKiller.hpp"

#include "PlaneBoundaryCondition.hpp"

#include "CryptShovingCaBasedDivisionRule.hpp"

#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "Warnings.hpp"

static const double M_END_TIME = 20.0; //2200

class TestCryptSimulationComparison : public AbstractCellBasedWithTimingsTestSuite
{
private:

    void GenerateStemCells(unsigned num_cells, std::vector<CellPtr>& rCells, double EquilibriumVolume)
    {
        double typical_cell_cycle_duration = 12.0;

        boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_cell_type(CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>());

        for (unsigned i=0; i<num_cells; i++)
        {
            ContactInhibitionCellCycleModel* p_model = new ContactInhibitionCellCycleModel();
            p_model->SetDimension(2);
            p_model->SetEquilibriumVolume(EquilibriumVolume);
            p_model->SetQuiescentVolumeFraction(0.8); //0.8 -> CI // 0.1 -> No CI!!!!

            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetCellProliferativeType(p_cell_type);
            double birth_time = - RandomNumberGenerator::Instance()->ranf() * typical_cell_cycle_duration;
            p_cell->SetBirthTime(birth_time);
            rCells.push_back(p_cell);
        }
    }

public:

    void NoTestVertexCrypt() throw (Exception)
    {
        double crypt_length = 12;
        double crypt_width = 5.5;

        // Create mesh
        unsigned cells_across = 6;

        CylindricalHoneycombVertexMeshGenerator generator(cells_across, 1, true);
        Cylindrical2dVertexMesh* p_mesh = generator.GetCylindricalMesh();

        p_mesh->Scale(crypt_width/cells_across,1.0);

        // Create cells
        std::vector<CellPtr> cells;
        GenerateStemCells(cells_across,cells,0.89); //1.0

        // Create tissue
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
        cell_population.AddCellWriter<CellVolumesWriter>();
        cell_population.AddCellWriter<CellIdWriter>();

        // Create crypt simulation from cell population
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetDt(1.0/200.0);
        simulator.SetSamplingTimestepMultiple(200);
        simulator.SetEndTime(M_END_TIME);
        simulator.SetOutputDirectory("CylindricalCrypt/Vertex");
        simulator.SetOutputDivisionLocations(true);
        simulator.SetOutputCellVelocities(true);

        // Add volume tracking modifier
        MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        // Create a force law and pass it to the simulation
        MAKE_PTR(NagaiHondaForce<2>, p_nagai_honda_force);
        simulator.AddForce(p_nagai_honda_force);

        // A NagaiHondaForce has to be used together with an AbstractTargetAreaModifier
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        // Solid base Boundary condition
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bcs, (&cell_population, zero_vector<double>(2), -unit_vector<double>(2,1)));
        p_bcs->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(p_bcs);

        // Sloughing killer
        MAKE_PTR_ARGS(PlaneBasedCellKiller<2>, p_killer, (&cell_population, crypt_length*unit_vector<double>(2,1), unit_vector<double>(2,1)));
        simulator.AddCellKiller(p_killer);

        // Run simulation
        simulator.Solve();

        Warnings::Instance()->QuietDestroy();
    }

    void NoTestMeshBasedCrypt() throw (Exception)
    {
        double crypt_length = 12-0.5;
        double crypt_width = 5.0;

        // Create mesh
        unsigned cells_across = 6;
        unsigned thickness_of_ghost_layer = 2;

        CylindricalHoneycombMeshGenerator generator(cells_across, 1, thickness_of_ghost_layer, crypt_width/cells_across);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        // Get location indices corresponding to real cells
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Create cells
        std::vector<CellPtr> cells;
        GenerateStemCells(cells_across,cells,0.8); //sqrt(3.0)/2.0

        // Create tissue
        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);
        cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
        cell_population.AddCellWriter<CellVolumesWriter>();
        cell_population.AddCellWriter<CellIdWriter>();

        // Create simulation from cell population
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetDt(1.0/200.0);
        simulator.SetSamplingTimestepMultiple(200);
        simulator.SetEndTime(M_END_TIME);
        simulator.SetOutputDirectory("CylindricalCrypt/Mesh");
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
        MAKE_PTR_ARGS(PlaneBasedCellKiller<2>, p_killer, (&cell_population, crypt_length*unit_vector<double>(2,1), unit_vector<double>(2,1)));
        simulator.AddCellKiller(p_killer);

        // Run simulation
        simulator.Solve();
    }

    void NoTestNodeBasedCrypt() throw (Exception)
    {
        double crypt_length = 12-0.5;
        double crypt_width = 5;

        // Create a simple mesh
        unsigned cells_across = 6;
        HoneycombMeshGenerator generator(cells_across, 1, 0, crypt_width/cells_across);
        TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();

        // Convert this to a Cylindrical2dNodesOnlyMesh
        Cylindrical2dNodesOnlyMesh* p_mesh = new Cylindrical2dNodesOnlyMesh(crypt_width);
        p_mesh->ConstructNodesWithoutMesh(*p_generating_mesh,crypt_width);

        // Get location indices corresponding to real cells
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Create cells
        std::vector<CellPtr> cells;
        GenerateStemCells(cells_across, cells, 0.825); // r=0.5 M_PI*0.25

        // Create a node-based cell population
        NodeBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
        cell_population.AddCellWriter<CellVolumesWriter>();
        cell_population.AddCellWriter<CellIdWriter>();

        for (unsigned index = 0; index < cell_population.rGetMesh().GetNumNodes(); index++)
        {
            cell_population.rGetMesh().GetNode(index)->SetRadius(0.5);
        }

        // Create simulation from cell population
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetDt(1.0/200.0);
        simulator.SetSamplingTimestepMultiple(200);
        simulator.SetEndTime(M_END_TIME);
        simulator.SetOutputDirectory("CylindricalCrypt/Node");
        simulator.SetOutputDivisionLocations(true);
        simulator.SetOutputCellVelocities(true);

        // Add volume tracking modifier
        MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        // Create a force law and pass it to the simulation
        MAKE_PTR(RepulsionForce<2>, p_repulsion_force);
        p_repulsion_force->SetMeinekeSpringStiffness(25.0);
        simulator.AddForce(p_repulsion_force);

        // Solid base boundary condition
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bcs, (&cell_population, zero_vector<double>(2), -unit_vector<double>(2,1)));
        p_bcs->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(p_bcs);

        // Sloughing killer
        MAKE_PTR_ARGS(PlaneBasedCellKiller<2>, p_killer, (&cell_population, crypt_length*unit_vector<double>(2,1), unit_vector<double>(2,1)));
        simulator.AddCellKiller(p_killer);

        // Run simulation
        simulator.Solve();

        // Clear memory
        delete p_mesh;
    }

    void NoTestPottsCrypt() throw (Exception)
    {
        double crypt_length = 12;
        double crypt_width = 6;
        unsigned cells_across = 6;

        // Create a simple 2D PottsMesh (periodic in x)
        PottsMeshGenerator<2> generator(crypt_width*4, cells_across, 4, (crypt_length+2)*4, 1, 4, 1, 1, 1, false, true);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        GenerateStemCells(cells_across,cells,16.0); //16.0

        // Create cell population
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetNumSweepsPerTimestep(10);
        cell_population.SetTemperature(0.001);
        cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
        cell_population.AddCellWriter<CellVolumesWriter>();
        cell_population.AddCellWriter<CellIdWriter>();

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("CylindricalCrypt/Potts");
        simulator.SetDt(0.1);
        simulator.SetSamplingTimestepMultiple(10);
        simulator.SetEndTime(M_END_TIME);
        simulator.SetOutputDivisionLocations(true);
        simulator.SetOutputCellVelocities(true);

        // Add volume tracking modifier
        MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        // Sloughing killer
        MAKE_PTR_ARGS(PlaneBasedCellKiller<2>, p_killer, (&cell_population, 4.0*crypt_length*unit_vector<double>(2,1), unit_vector<double>(2,1)));
        simulator.AddCellKiller(p_killer);

        // Create update rules and pass to the simulation
        MAKE_PTR(VolumeConstraintPottsUpdateRule<2>, p_volume_constraint_update_rule);
        simulator.AddPottsUpdateRule(p_volume_constraint_update_rule);
        MAKE_PTR(AdhesionPottsUpdateRule<2>, p_adhesion_update_rule);
        simulator.AddPottsUpdateRule(p_adhesion_update_rule);

        // Run simulation
        simulator.Solve();
    }

    void TestCaCrypt() throw (Exception)
    {
        unsigned cells_across = 6;
        unsigned crypt_length = 12;

        // Create a simple 2D PottsMesh (periodic in x)
        PottsMeshGenerator<2> generator(cells_across, 0, 0, crypt_length*2+1, 0, 0, 1, 0, 0, false, true);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Specify where cells lie
        std::vector<unsigned> location_indices;
        for (unsigned index=0; index<cells_across; index++)
        {
            location_indices.push_back(index);
        }


//        // Create cells
//        std::vector<CellPtr> cells;
//        GenerateStemCells(cells_across,cells,1.0);

		MAKE_PTR(TransitCellProliferativeType, p_transit_type);

		std::vector<CellPtr> cells;
		CellsGenerator<StochasticDurationCellCycleModel, 2> cells_generator;
		cells_generator.GenerateBasicRandom(cells, location_indices.size(), p_transit_type);


        // Create cell population
        CaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);
        cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
        cell_population.AddCellWriter<CellVolumesWriter>();
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellAncestorWriter>();

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("CylindricalCrypt/Ca");
        simulator.SetDt(1);
        simulator.SetEndTime(100);//M_END_TIME);
        simulator.SetOutputDivisionLocations(true);
        simulator.SetOutputCellVelocities(true);

        // Add Volume tracking modifier
        MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        // Add Division Rule
        boost::shared_ptr<AbstractCaBasedDivisionRule<2> > p_division_rule(new CryptShovingCaBasedDivisionRule());
        cell_population.SetCaBasedDivisionRule(p_division_rule);



        // Sloughing killer
        MAKE_PTR_ARGS(PlaneBasedCellKiller<2>, p_killer, (&cell_population, crypt_length*unit_vector<double>(2,1), unit_vector<double>(2,1)));
        simulator.AddCellKiller(p_killer);

        // Run simulation
        simulator.Solve();

        simulator.SetEndTime(500);

        simulator.rGetCellPopulation().SetCellAncestorsToLocationIndices();

        // Run simulation to new end time
        simulator.Solve();



    }
};

#endif /*TESTCRYPTSIMULATIONCOMPARISON_HPP_*/
