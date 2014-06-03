
#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "AbstractCellBasedTestSuite.hpp"
#include "CellLabel.hpp"
#include "SmartPointers.hpp"
#include "CellsGenerator.hpp"
#include "DeltaNotchCellCycleModel.hpp"
#include "WildTypeCellMutationState.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

#include "OffLatticeSimulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "NagaiHondaDifferentialAdhesionForce.hpp"
#include "VertexDiffusionForce.hpp"
#include "SimpleTargetAreaModifier.hpp"

#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "DifferentialAdhesionSpringForce.hpp"
#include "SensibleDiffusionForce.hpp"
#include "RepulsionForce.hpp"
#include "GeneralisedLinearSpringForce.hpp"

#include "OnLatticeSimulation.hpp"
#include "DeltaNotchTrackingModifier.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "PottsMeshGenerator.hpp"
#include "VolumeConstraintPottsUpdateRule.hpp"
#include "AdhesionPottsUpdateRule.hpp"
#include "DifferentialAdhesionPottsUpdateRule.hpp"

#include "CellIdWriter.hpp"
#include "DeltaNotchWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "VoronoiDataWriter.hpp"

#include "PetscSetupAndFinalize.hpp"

static const double M_TIME_FOR_SIMULATION = 100.0;

static const double M_NUM_CELLS_ACROSS = 10; // this ^2 cells



class TestDeltaNotch: public AbstractCellBasedTestSuite
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

    void GenerateCells(unsigned num_cells, std::vector<CellPtr>& rCells)
    {
    	double typical_transit_cell_cycle_duration = 12.0;

        boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
        boost::shared_ptr<AbstractCellProperty> p_prolif_type(CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>());
//        boost::shared_ptr<AbstractCellProperty> p_prolif_type(CellPropertyRegistry::Instance()->Get<TransitCellProliferativeType>());

        for (unsigned i=0; i<num_cells; i++)
        {
            std::vector<double> initial_conditions;
            initial_conditions.push_back(RandomNumberGenerator::Instance()->ranf());
            initial_conditions.push_back(RandomNumberGenerator::Instance()->ranf());

            DeltaNotchCellCycleModel* p_model = new DeltaNotchCellCycleModel();
            p_model->SetInitialConditions(initial_conditions);
            p_model->SetDimension(2);
            p_model->SetMaxTransitGenerations(2u);
            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetCellProliferativeType(p_prolif_type);
			double birth_time = SimulationTime::Instance()->GetTime() - RandomNumberGenerator::Instance()->ranf() * typical_transit_cell_cycle_duration;
			p_cell->SetBirthTime(birth_time);
            rCells.push_back(p_cell);
        }
    }

public:

    void TestVertexMonolayerDeltaNotch() throw (Exception)
    {
        // Create a simple 2D MutableVertexMesh
        HoneycombVertexMeshGenerator generator(M_NUM_CELLS_ACROSS,M_NUM_CELLS_ACROSS);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();
        p_mesh->SetCellRearrangementThreshold(0.1);

        // Slows things down but can use a larger timestep and diffusion forces.
        p_mesh->SetCheckForInternalIntersections(true);

        // ASsociate each cell with a cell-cycle model that incorporates a Delta-Notch ODE system
        std::vector<CellPtr> cells;
        GenerateCells(p_mesh->GetNumElements(),cells);

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set population to output all data to results files
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<DeltaNotchWriter>();

        // Set up cell-based simulation and output directory
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("DeltaNotch/Vertex");

        // Set time step and end time for simulation
        simulator.SetDt(0.1);
		simulator.SetSamplingTimestepMultiple(10);
        simulator.SetEndTime(M_TIME_FOR_SIMULATION);

        // Add DeltaNotch modifier
        MAKE_PTR(DeltaNotchTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        // Create Forces and pass to simulation NOTE: PARAMETERS CHOSEN TO GET CIRCULAR MONOLAYER
        MAKE_PTR(NagaiHondaForce<2>, p_force);
        p_force->SetNagaiHondaDeformationEnergyParameter(5.5);
		p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(0.0);
		p_force->SetNagaiHondaCellCellAdhesionEnergyParameter(0.6);
		p_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(1.2);
        simulator.AddForce(p_force);

        // A NagaiHondaForce has to be used together with an AbstractTargetAreaModifier
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);
        // Run simulation
        simulator.Solve();

        // Check that the same number of cells
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), 100u);

        // Test no births or deaths
        TS_ASSERT_EQUALS(simulator.GetNumBirths(), 0u);
        TS_ASSERT_EQUALS(simulator.GetNumDeaths(), 0u);

   }


    void TestPottsMonolayerDeltaNotch() throw (Exception)
    {
        // Create a simple 2D PottsMesh
        unsigned element_size = 4;
        unsigned domain_size = 158;//M_NUM_CELLS_ACROSS * element_size * 2; // Twice the initial domain size
        PottsMeshGenerator<2> generator(domain_size, M_NUM_CELLS_ACROSS, element_size, domain_size, M_NUM_CELLS_ACROSS, element_size);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        GenerateCells(p_mesh->GetNumElements(),cells);

        // Create cell population
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<DeltaNotchWriter>();
        //cell_population.SetNumSweepsPerTimestep(10);
       // cell_population.SetTemperature(0.01);// This is the default value


        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("DeltaNotch/Potts");

        // Set time step and end time for simulation
        simulator.SetDt(0.1); // This is the default value
        simulator.SetSamplingTimestepMultiple(10);
        simulator.SetEndTime(M_TIME_FOR_SIMULATION);

        // Add DeltaNotch modifier
        MAKE_PTR(DeltaNotchTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        // Create update rules and pass to the simulation
        MAKE_PTR(VolumeConstraintPottsUpdateRule<2>, p_volume_constraint_update_rule);
        p_volume_constraint_update_rule->SetMatureCellTargetVolume(16); // i.e 4x4 cells
        p_volume_constraint_update_rule->SetDeformationEnergyParameter(0.2);
        simulator.AddPottsUpdateRule(p_volume_constraint_update_rule);

        MAKE_PTR(AdhesionPottsUpdateRule<2>, p_adhesion_update_rule);
        simulator.AddPottsUpdateRule(p_adhesion_update_rule);

        // Run simulation
        simulator.Solve();

        // Check that the same number of cells
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), 100u);

        // Test no births or deaths
        TS_ASSERT_EQUALS(simulator.GetNumBirths(), 0u);
        TS_ASSERT_EQUALS(simulator.GetNumDeaths(), 0u);
    }


    void TestMeshBasedWithGhostsMonolayerDeltaNotch() throw (Exception)
	{
		// Create a simple mesh
		unsigned num_ghosts = 3;
		HoneycombMeshGenerator generator(M_NUM_CELLS_ACROSS, M_NUM_CELLS_ACROSS, num_ghosts);
		MutableMesh<2,2>* p_mesh = generator.GetMesh();

		// Set up cells, one for each non ghost Node
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();//**Changed**//
        std::vector<CellPtr> cells;
        GenerateCells(location_indices.size(),cells);

		// Create cell population
		MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);

		// Set population to output all data to results files
		cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<DeltaNotchWriter>();
        cell_population.SetWriteVtkAsPoints(false);
        cell_population.AddPopulationWriter<VoronoiDataWriter>();

		// Set up cell-based simulation and output directory
		OffLatticeSimulation<2> simulator(cell_population);
		simulator.SetOutputDirectory("DeltaNotch/MeshGhost");

		// Set time step and end time for simulation
		simulator.SetDt(0.01);
		simulator.SetSamplingTimestepMultiple(100);
		simulator.SetEndTime(M_TIME_FOR_SIMULATION);

        // Add DeltaNotch modifier
        MAKE_PTR(DeltaNotchTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        simulator.AddForce(p_force);

      	// Run simulation
		simulator.Solve();

		// Check that the same number of cells
		TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), 100u);

		// Test no births or deaths
		TS_ASSERT_EQUALS(simulator.GetNumBirths(), 0u);
		TS_ASSERT_EQUALS(simulator.GetNumDeaths(), 0u);
   }

    void TestMeshBasedMonolayerDeltaNotch() throw (Exception)
	{
        // Create a simple mesh
        unsigned num_ghosts = 0;
        HoneycombMeshGenerator generator(M_NUM_CELLS_ACROSS, M_NUM_CELLS_ACROSS, num_ghosts);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

		// Set up cells, one for each Node
        std::vector<CellPtr> cells;
        GenerateCells(p_mesh->GetNumNodes(),cells);

		// Create cell population
		MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

		// Set population to output all data to results files
		cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<DeltaNotchWriter>();
        cell_population.SetWriteVtkAsPoints(false);
		cell_population.AddPopulationWriter<VoronoiDataWriter>();

		// Set up cell-based simulation and output directory
		OffLatticeSimulation<2> simulator(cell_population);
		simulator.SetOutputDirectory("DeltaNotch/Mesh");

		// Set time step and end time for simulation
		simulator.SetDt(0.01);
		simulator.SetSamplingTimestepMultiple(100);
		simulator.SetEndTime(M_TIME_FOR_SIMULATION);

        // Add DeltaNotch modifier
        MAKE_PTR(DeltaNotchTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);


        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        simulator.AddForce(p_force);

      	// Run simulation
		simulator.Solve();

		// Check that the same number of cells
		TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), 100u);

		// Test no births or deaths
		TS_ASSERT_EQUALS(simulator.GetNumBirths(), 0u);
		TS_ASSERT_EQUALS(simulator.GetNumDeaths(), 0u);
    }


    void TestNodeBasedMonolayerCellSorting() throw (Exception)
	{
        // Create a simple mesh
        HoneycombMeshGenerator generator(M_NUM_CELLS_ACROSS, M_NUM_CELLS_ACROSS, 0);
        TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);

		// Set up cells, one for each Node
		std::vector<CellPtr> cells;
        GenerateCells(mesh.GetNumNodes(),cells);

		// Create cell population
		NodeBasedCellPopulation<2> cell_population(mesh, cells);

		// Set population to output all data to results files
		cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<DeltaNotchWriter>();
		cell_population.AddPopulationWriter<VoronoiDataWriter>();

		// Set up cell-based simulation and output directory
		OffLatticeSimulation<2> simulator(cell_population);
		simulator.SetOutputDirectory("DeltaNotch/Node");

		// Set time step and end time for simulation
		simulator.SetDt(0.01);
		simulator.SetSamplingTimestepMultiple(100);
		simulator.SetEndTime(M_TIME_FOR_SIMULATION);

        // Add DeltaNotch modifier
        MAKE_PTR(DeltaNotchTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        MAKE_PTR(RepulsionForce<2>, p_force);
        simulator.AddForce(p_force);

		// Run simulation
		simulator.Solve();

		// Check that the same number of cells
		TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), 100u);

		// Test no births or deaths
		TS_ASSERT_EQUALS(simulator.GetNumBirths(), 0u);
		TS_ASSERT_EQUALS(simulator.GetNumDeaths(), 0u);
   }


};
