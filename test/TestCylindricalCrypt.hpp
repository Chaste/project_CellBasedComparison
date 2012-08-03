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

#include "PottsMeshGenerator.hpp"
#include "SimpleWntCellCycleModel.hpp"
#include "OnLatticeSimulation.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "VolumeConstraintPottsUpdateRule.hpp"
#include "AdhesionPottsUpdateRule.hpp"

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

    void TestSimpleVertexCryptforComparison() throw (Exception)
    {
        double crypt_length = 12;

        // Create mesh
        unsigned cells_across = 6;
        unsigned cells_up = 14;
        CylindricalHoneycombVertexMeshGenerator generator(cells_across, cells_up);
        Cylindrical2dVertexMesh* p_mesh = generator.GetCylindricalMesh();

        // Set up cells
        std::vector<CellPtr> cells;
        CryptCellsGenerator<StochasticDurationGenerationBasedCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, std::vector<unsigned>(), true, 0.5, 2.0, 4.0, 5.8);

        // Create tissue
        VertexBasedCellPopulation<2> crypt(*p_mesh, cells);

        // Create crypt simulation from cell population
        OffLatticeSimulation<2> simulator(crypt);
        simulator.SetSamplingTimestepMultiple(50);
        simulator.SetEndTime(1);
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

    void TestSimpleMeshBasedCryptforComparison() throw (Exception)
    {
        double crypt_length = 12;

        // Create mesh
        unsigned cells_across = 6;
        unsigned cells_up = 12;
        double crypt_width = 5.0;
        unsigned thickness_of_ghost_layer = 1;

        CylindricalHoneycombMeshGenerator generator(cells_across, cells_up, thickness_of_ghost_layer, crypt_width/cells_across);
        Cylindrical2dMesh* p_mesh = generator.GetCylindricalMesh();

        // Get location indices corresponding to real cells
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Set up cells
        std::vector<CellPtr> cells;
        CryptCellsGenerator<StochasticDurationGenerationBasedCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, location_indices, true, 0.5, 2.0, 4.0, 5.8);

        // Create tissue
        MeshBasedCellPopulationWithGhostNodes<2> crypt(*p_mesh, cells, location_indices);

        // Create crypt simulation from cell population
        OffLatticeSimulation<2> simulator(crypt);
        simulator.SetSamplingTimestepMultiple(12);
        simulator.SetEndTime(10);
        simulator.SetOutputDirectory("CylindricalCrypt/Mesh");

        // Add Crypt BCS
        MAKE_PTR_ARGS(CryptSimulationBoundaryCondition<2>, p_crypt_bcs,(&crypt));
        p_crypt_bcs->SetUseJiggledBottomCells(true);
        simulator.AddCellPopulationBoundaryCondition(p_crypt_bcs);

        // Create a force law and pass it to the simulation
        GeneralisedLinearSpringForce<2> linear_force;
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        simulator.AddForce(p_linear_force);

        // Make crypt shorter for sloughing
        MAKE_PTR_ARGS(SloughingCellKiller<2>, p_killer, (&crypt, crypt_length));
        simulator.AddCellKiller(p_killer);

        // Run simulation
        TS_ASSERT_THROWS_NOTHING(simulator.Solve());
    }


    void TestNodeBasedCrypt() throw (Exception)
    {
        double crypt_length = 12;

        // Create a simple mesh
        unsigned num_cells_depth = 14;
        unsigned num_cells_width = 5;
        HoneycombMeshGenerator generator(num_cells_width, num_cells_depth, 0);
        TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh);

        // Get location indices corresponding to real cells
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Create cells
        std::vector<CellPtr> cells;
        CryptCellsGenerator<StochasticDurationGenerationBasedCellCycleModel> cells_generator;
        cells_generator.Generate(cells, &mesh, location_indices, true, 0.5, 2.0, 4.0, 5.8);


        for (unsigned index = 0; index < cells.size(); index++)
        {
            //dynamic_cast<AbstractSimpleGenerationBasedCellCycleModel*>(cells[index]->GetCellCycleModel())->SetMaxTransitGenerations(2u);
        }

        // Create a node-based cell population
        NodeBasedCellPopulation<2> crypt(mesh, cells);
        crypt.SetMechanicsCutOffLength(1.5);

        // Create crypt simulation from cell population
        OffLatticeSimulation<2> simulator(crypt);
        simulator.SetDt(0.001);
        simulator.SetSamplingTimestepMultiple(100);
        simulator.SetEndTime(1.0);
        simulator.SetOutputDirectory("CylindricalCrypt/Node");

        // Create a force law and pass it to the simulation
        GeneralisedLinearSpringForce<2> linear_force;
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        simulator.AddForce(p_linear_force);

        // Make crypt shorter for sloughing
        MAKE_PTR_ARGS(SloughingCellKiller<2>, p_killer, (&crypt, crypt_length));
        simulator.AddCellKiller(p_killer);

        // Run simulation
        simulator.Solve();
    }

    void TestPottsCrypt() throw (Exception)
    {
        double crypt_length = 40;

        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(20, 5, 4, 45, 10, 4, 1, 1, 1, true);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        CryptCellsGenerator<StochasticDurationGenerationBasedCellCycleModel> cells_generator;
        cells_generator.Generate(cells, p_mesh, std::vector<unsigned>(), true, 2.5, 8.0, 16.0, 36);

        // Create cell population
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetOutputCellVolumes(true);

        // Create an instance of a Wnt concentration
//        WntConcentration<2>::Instance()->SetType(LINEAR);
//        WntConcentration<2>::Instance()->SetCellPopulation(cell_population);
//        WntConcentration<2>::Instance()->SetCryptLength(crypt_length);

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("CylindricalCrypt/Potts");
        simulator.SetDt(0.1);
        simulator.SetSamplingTimestepMultiple(1);
        simulator.SetEndTime(10.0);
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
