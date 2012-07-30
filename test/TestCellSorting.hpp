
#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "AbstractCellBasedTestSuite.hpp"
#include "CellLabel.hpp"
#include "SmartPointers.hpp"
#include "CellsGenerator.hpp"
#include "StochasticDurationCellCycleModel.hpp"

#include "OffLatticeSimulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "NagaiHondaDifferentialAdhesionForce.hpp"
#include "VertexDiffusionForce.hpp"

#include "OnLatticeSimulation.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "PottsMeshGenerator.hpp"
#include "VolumeConstraintPottsUpdateRule.hpp"
#include "AdhesionPottsUpdateRule.hpp"
#include "DifferentialAdhesionPottsUpdateRule.hpp"

class TestCellSorting: public AbstractCellBasedTestSuite
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

    void RandomlyLabelCells(std::vector<CellPtr>& rCells, boost::shared_ptr<AbstractCellProperty> pLabel, double labelledRatio)
    {
        for (unsigned i = 0; i<rCells.size(); i++)
        {
            if (RandomNumberGenerator::Instance()->ranf() < labelledRatio)
            {
                rCells[i]->AddCellProperty(pLabel);
            }
        }
    }

public:

    /**
     * Simulate a population of cells exhibiting cell sorting using the
     * vertex dynamics model proposed by T. Nagai and H. Honda ("A dynamic
     * cell model for the formation of epithelial tissues", Philosophical
     * Magazine Part B 81:699-719).
     *
     * Each of the vertex dynamics model parameter member variables are
     * rescaled such that mDampingConstantNormal takes the default value 1,
     * whereas Nagai and Honda (who denote the parameter by nu) take the
     * value 0.01.
     */
    void TestVertexMonolayerCellSorting() throw (Exception)
    {
        // Create a simple 2D MutableVertexMesh
        HoneycombVertexMeshGenerator generator(10, 10);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();
        p_mesh->SetCellRearrangementThreshold(0.1);

        // Set up cells, one for each VertexElement
        std::vector<CellPtr> cells;
        CellsGenerator<StochasticDurationCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells,p_mesh->GetNumElements(), DIFFERENTIATED);

        // Randomly label some cells
        boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<CellLabel>());
        RandomlyLabelCells(cells, p_state, 0.5);

        // Create cell population
        VertexBasedCellPopulation<2> population(*p_mesh, cells);

        // Set population to output all data to results files
        population.SetOutputCellIdData(true);
        population.SetOutputCellMutationStates(true);
        population.SetOutputCellAncestors(true);
        population.SetOutputCellProliferativeTypes(true);
        population.SetOutputCellVariables(true);
        population.SetOutputCellCyclePhases(true);
        population.SetOutputCellAges(true);
        population.SetOutputCellVolumes(true);

        // Set up cell-based simulation and output directory
        OffLatticeSimulation<2> simulator(population);
        simulator.SetOutputDirectory("CellSorting/Vertex");

        // Set time step and end time for simulation
        simulator.SetDt(0.001);
        simulator.SetEndTime(10.0);

        // Only record results every 100 time steps
        simulator.SetSamplingTimestepMultiple(100);

        // Set up force law and pass it to the simulation
        MAKE_PTR(NagaiHondaDifferentialAdhesionForce<2>, p_force);
        p_force->SetNagaiHondaDeformationEnergyParameter(55.0);
        p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(0.0);
        p_force->SetNagaiHondaCellCellAdhesionEnergyParameter(1.0);
        p_force->SetNagaiHondaLabeledCellCellAdhesionEnergyParameter(6.0);
        p_force->SetNagaiHondaLabeledCellLabeledCellAdhesionEnergyParameter(3.0);
        p_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(12.0);
        p_force->SetNagaiHondaLabeledCellBoundaryAdhesionEnergyParameter(40.0);
        simulator.AddForce(p_force);

//      you need to check for internal intersections if running with diffusion.
//        MAKE_PTR_ARGS(VertexDiffusionForce<2>, p_random_force, (0.01));
//        simulator.AddForce(p_random_force);

        // Run simulation
        simulator.Solve();

        // Check that the same number of cells
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), 100u);

        // Test no births or deaths
        TS_ASSERT_EQUALS(simulator.GetNumBirths(), 0u);
        TS_ASSERT_EQUALS(simulator.GetNumDeaths(), 0u);

   }


    /**
     * Simulate a population of cells exhibiting cell sorting using the
     * Potts model.
     */
    void TestPottsMonolayerCellSorting() throw (Exception)
    {
        // Create a simple 3D PottsMesh
        unsigned domain_size = 80;
        unsigned element_number = 10;
        unsigned element_size = 4;
        PottsMeshGenerator<2> generator(domain_size, element_number, element_size, domain_size, element_number, element_size);
        PottsMesh<2>* p_mesh = generator.GetMesh();


        // Create cells
        std::vector<CellPtr> cells;
        CellsGenerator<StochasticDurationCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), DIFFERENTIATED);

        // Randomly label some cells
        boost::shared_ptr<AbstractCellProperty> p_label(CellPropertyRegistry::Instance()->Get<CellLabel>());
        RandomlyLabelCells(cells, p_label, 0.5);

        // Create cell population
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetOutputCellMutationStates(true); // So outputs the labelled cells
        cell_population.SetNumSweepsPerTimestep(1); // This is the default value

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("CellSorting/Potts");
        simulator.SetDt(0.1); // This is the default value
        simulator.SetEndTime(100.0); // i.e 1000 MCS

        // Create update rules and pass to the simulation
        MAKE_PTR(VolumeConstraintPottsUpdateRule<2>, p_volume_constraint_update_rule);
        p_volume_constraint_update_rule->SetMatureCellTargetVolume(16); // i.e 4x4 cells
        p_volume_constraint_update_rule->SetDeformationEnergyParameter(0.2);
        simulator.AddPottsUpdateRule(p_volume_constraint_update_rule);

        MAKE_PTR(DifferentialAdhesionPottsUpdateRule<2>, p_differential_adhesion_update_rule);
        p_differential_adhesion_update_rule->SetLabelledCellLabelledCellAdhesionEnergyParameter(0.16);
        p_differential_adhesion_update_rule->SetLabelledCellCellAdhesionEnergyParameter(0.11);
        p_differential_adhesion_update_rule->SetCellCellAdhesionEnergyParameter(0.02);
        p_differential_adhesion_update_rule->SetLabelledCellBoundaryAdhesionEnergyParameter(0.16);
        p_differential_adhesion_update_rule->SetCellBoundaryAdhesionEnergyParameter(0.16);
        simulator.AddPottsUpdateRule(p_differential_adhesion_update_rule);

        // Run simulation
        simulator.Solve();

        // Check that the same number of cells
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), 100u);

        // Test no births or deaths
        TS_ASSERT_EQUALS(simulator.GetNumBirths(), 0u);
        TS_ASSERT_EQUALS(simulator.GetNumDeaths(), 0u);
    }


};
