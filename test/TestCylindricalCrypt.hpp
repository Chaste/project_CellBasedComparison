/*
 * Would prefer to use Off and On Lattice sim for all these then can discuss added features of Crypt Sim 2D
 *
 */
#ifndef TESTCRYPTSIMULATIONCOMPARISON_HPP_
#define TESTCRYPTSIMULATIONCOMPARISON_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before any other cell_based headers
#include "CellBasedSimulationArchiver.hpp"
#include "CryptSimulation2d.hpp"
//#include "NodeBasedCryptSimulation2d.hpp"
#include "CryptCellsGenerator.hpp"
#include "CellsGenerator.hpp"
//#include "VanLeeuwen2009WntSwatCellCycleModelHypothesisOne.hpp"
#include "LinearSpringWithVariableSpringConstantsForce.hpp"
#include "CylindricalHoneycombVertexMeshGenerator.hpp"
#include "CylindricalHoneycombMeshGenerator.hpp"
//#include "TargetedCellKiller.hpp"
#include "RandomCellKiller.hpp"
#include "SloughingCellKiller.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "NumericFileComparison.hpp"
#include "CellBasedEventHandler.hpp"
#include "RepulsionForce.hpp"
#include "StochasticDurationGenerationBasedCellCycleModel.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "NagaiHondaForce.hpp"
#include "Cylindrical2dNodesOnlyMesh.hpp"

#include "PeriodicNodeBasedBoundaryCondition.hpp"
#include "StemCellRetainerForce.hpp"

#include "PottsMeshGenerator.hpp"
#include "SimpleWntCellCycleModel.hpp"
#include "OnLatticeSimulation.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "VolumeConstraintPottsUpdateRule.hpp"
#include "AdhesionPottsUpdateRule.hpp"

#include "Debug.hpp"

class TestCryptSimulationcomparison : public AbstractCellBasedTestSuite
{
private:


    double mLastStartTime;
    void setUp()
    {
        mLastStartTime = std::clock();
        AbstractCellBasedTestSuite::setUp();
    }
    void tearDown()
    {
        double time = std::clock();
        double elapsed_time = (time - mLastStartTime)/(CLOCKS_PER_SEC);
        std::cout << "Elapsed time: " << elapsed_time << std::endl;
        AbstractCellBasedTestSuite::tearDown();
    }

public:

    void nTestSimpleVertexCryptforComparison() throw (Exception)
    {
        double crypt_length = 12;
        double crypt_width = 6.0;

        // Create mesh
        unsigned cells_across = 6;
        unsigned cells_up = 14;
        CylindricalHoneycombVertexMeshGenerator generator(cells_across, cells_up, crypt_width/cells_across);
        Cylindrical2dVertexMesh* p_mesh = generator.GetCylindricalMesh();

        // Set up cells
        std::vector<CellPtr> cells;
        CryptCellsGenerator<StochasticDurationGenerationBasedCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, std::vector<unsigned>(), true, 1.0, 2.5, 4.5, 6.0);

        // Create tissue
        VertexBasedCellPopulation<2> crypt(*p_mesh, cells);

        // Create crypt simulation from cell population
        OffLatticeSimulation<2> simulator(crypt);
        simulator.SetSamplingTimestepMultiple(50);
        simulator.SetEndTime(100);
        simulator.SetOutputDirectory("CylindricalCrypt/Vertex");

        // Add Crypt BCS
        MAKE_PTR_ARGS(CryptSimulationBoundaryCondition<2>, p_crypt_bcs,(&crypt));
        p_crypt_bcs->SetUseJiggledBottomCells(true);
        simulator.AddCellPopulationBoundaryCondition(p_crypt_bcs);

        // Create a force law and pass it to the simulation
        MAKE_PTR(NagaiHondaForce<2>, p_nagai_honda_force);
        simulator.AddForce(p_nagai_honda_force);

        // Make crypt shorter for sloughing
        MAKE_PTR_ARGS(SloughingCellKiller<2>, p_killer, (&crypt, crypt_length));
        simulator.AddCellKiller(p_killer);

        // Run simulation
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());
    }

    void nTestSimpleMeshBasedCryptforComparison() throw (Exception)
    {
        double crypt_length = 12-0.5;
        double crypt_width = 6.0;

        // Create mesh
        unsigned cells_across = 6;
        unsigned cells_up = 14;
        unsigned thickness_of_ghost_layer = 1;

        CylindricalHoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer, crypt_width/cells_across);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        // Get location indices corresponding to real cells
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Set up cells
        std::vector<CellPtr> cells;
        CryptCellsGenerator<StochasticDurationGenerationBasedCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, location_indices, true, 1.0-0.5, 2.5-0.5, 4.5- 0.5, 6.0 -0.5);

        // Create tissue
        MeshBasedCellPopulationWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);

        // Create crypt simulation from cell population
        OffLatticeSimulation<2> simulator(crypt);
        simulator.SetSamplingTimestepMultiple(12);
        simulator.SetEndTime(100);
        simulator.SetOutputDirectory("CylindricalCrypt/Mesh");

        // Add Crypt BCS
        MAKE_PTR_ARGS(CryptSimulationBoundaryCondition<2>, p_crypt_bcs,(&crypt));
        p_crypt_bcs->SetUseJiggledBottomCells(true);
        simulator.AddCellPopulationBoundaryCondition(p_crypt_bcs);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetMeinekeSpringStiffness(30.0); //normally 15.0 but 30 in all CellBased Papers Modified to stop crowding at base of crypt        simulator.AddForce(p_linear_force);
        simulator.AddForce(p_linear_force);

        // Make crypt shorter for sloughing
        MAKE_PTR_ARGS(SloughingCellKiller<2>, p_killer, (&crypt, crypt_length));
        simulator.AddCellKiller(p_killer);

        // Run simulation
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());
    }


    void TestNodeBasedCrypt() throw (Exception)
    {
        double crypt_length = 12-0.5;
        double crypt_width = 6.0;

        // Create a simple mesh
        unsigned cells_across = 6;
        unsigned cells_up = 14;
        HoneycombMeshGenerator generator(cells_across, cells_up, 0, crypt_width/cells_across);
        TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();

        // Convert this to a Cylindrical2dNodesOnlyMesh
        Cylindrical2dNodesOnlyMesh* p_mesh = new Cylindrical2dNodesOnlyMesh(crypt_width);
        p_mesh->ConstructNodesWithoutMesh(*p_generating_mesh);

        // Get location indices corresponding to real cells
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Create cells
        std::vector<CellPtr> cells;
        CryptCellsGenerator<StochasticDurationGenerationBasedCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, location_indices, true, 1.0-0.5, 2.5-0.5, 4.5- 0.5, 6.0 -0.5);
        //cells_generator.Generate(cells, p_mesh, location_indices, true, 0.0, 0.0, 0.0,0.0,0.0);

        // Create a node-based cell population
        NodeBasedCellPopulation<2> crypt(*p_mesh, cells);
        crypt.SetMechanicsCutOffLength(crypt_width/2.0);

        for (unsigned index = 0; index < crypt.rGetMesh().GetNumNodes(); index++)
        {
            PRINT_VARIABLE(index);
            crypt.rGetMesh().SetCellRadius(index,0.53);
        }
        TS_ASSERT_DELTA(crypt.rGetMesh().GetCellRadius(0),0.53,1e-6);

        // Create crypt simulation from cell population
        OffLatticeSimulation<2> simulator(crypt);
        //simulator.SetDt(0.001);
        simulator.SetSamplingTimestepMultiple(12);
        simulator.SetEndTime(100.0);
        simulator.SetOutputDirectory("CylindricalCrypt/Node");

        // Create a force law and pass it to the simulation
        MAKE_PTR(RepulsionForce<2>, p_repulsion_force);
        p_repulsion_force->SetMeinekeSpringStiffness(30.0); //normally 15.0 but 30 in all CellBased Papers Modified to stop crowding at base of crypt
        simulator.AddForce(p_repulsion_force);

        // Make crypt shorter for sloughing
        MAKE_PTR_ARGS(SloughingCellKiller<2>, p_killer, (&crypt, crypt_length));
        simulator.AddCellKiller(p_killer);

        // Add Periodic BCS to move nodes round periodic boundary
        MAKE_PTR_ARGS(PeriodicNodeBasedBoundaryCondition, p_bcs,(&crypt, crypt_width));
        simulator.AddCellPopulationBoundaryCondition(p_bcs);

        // Add Crypt BCS
        MAKE_PTR_ARGS(CryptSimulationBoundaryCondition<2>, p_crypt_bcs,(&crypt));
        p_crypt_bcs->SetUseJiggledBottomCells(true);
        simulator.AddCellPopulationBoundaryCondition(p_crypt_bcs);

        // Create a force law to retain stem cells niche and pass it to the simulation
        MAKE_PTR(StemCellRetainerForce<2>, p_retainer_force);
        p_retainer_force->SetForceMagnitudeParameter(50.0);
        simulator.AddForce(p_retainer_force);

        // Run simulation
        simulator.Solve();

        // Clear memory
        delete p_mesh;
    }

    void TestPottsCrypt() throw (Exception)
    {
        double crypt_length = 12*4;

        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(6*4, 6, 4, 13*4, 12, 4, 1, 1, 1, true);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CryptCellsGenerator<StochasticDurationGenerationBasedCellCycleModel> cells_generator;
        //cells_generator.Generate(cells, p_mesh, std::vector<unsigned>(), true, 2.5, 8.0, 16.0, 36);
        cells_generator.Generate(cells, p_mesh, std::vector<unsigned>(), true, 4.0*1.0, 4.0*2.5, 4.0*4.5, 4.0*6.0);

        // Create cell population
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetOutputCellVolumes(true);
        cell_population.SetNumSweepsPerTimestep(10);
        cell_population.SetTemperature(0.01);

        // Create an instance of a Wnt concentration
//        WntConcentration<2>::Instance()->SetType(LINEAR);
//        WntConcentration<2>::Instance()->SetCellPopulation(cell_population);
//        WntConcentration<2>::Instance()->SetCryptLength(crypt_length);

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("CylindricalCrypt/Potts");
        simulator.SetDt(0.1);
        simulator.SetSamplingTimestepMultiple(1);
        simulator.SetEndTime(100.0);
        simulator.SetOutputCellVelocities(true);

        // Create cell killer and pass in to simulation
        MAKE_PTR_ARGS(SloughingCellKiller<2>, p_killer, (&cell_population, crypt_length));
        simulator.AddCellKiller(p_killer);

        // Create update rules and pass to the simulation
        MAKE_PTR(VolumeConstraintPottsUpdateRule<2>, p_volume_constraint_update_rule);
        simulator.AddPottsUpdateRule(p_volume_constraint_update_rule);
        MAKE_PTR(AdhesionPottsUpdateRule<2>, p_adhesion_update_rule);
        simulator.AddPottsUpdateRule(p_adhesion_update_rule);

        // Run simulation
        simulator.Solve();
    }
};

#endif /*TESTCRYPTSIMULATIONCOMPARISON_HPP_*/
