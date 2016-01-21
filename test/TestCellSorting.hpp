
#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "CellLabel.hpp"
#include "SmartPointers.hpp"
#include "CellsGenerator.hpp"
#include "StochasticDurationGenerationBasedCellCycleModel.hpp"
#include "TransitCellProliferativeType.hpp"

#include "HeterotypicBoundaryLengthWriter.hpp"

#include "OffLatticeSimulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "NagaiHondaDifferentialAdhesionForce.hpp"
#include "RandomMotionForce.hpp"
#include "SimpleTargetAreaModifier.hpp"

#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "DifferentialAdhesionGeneralisedLinearSpringForce.hpp"
#include "RandomMotionForce.hpp"

#include "OnLatticeSimulation.hpp"
#include "CellPopulationAdjacencyMatrixWriter.hpp"

#include "CaBasedCellPopulation.hpp"
#include "ShovingCaBasedDivisionRule.hpp"
#include "RandomCaSwitchingUpdateRule.hpp"
#include "DifferentialAdhesionCaSwitchingUpdateRule.hpp"

#include "PottsBasedCellPopulation.hpp"
#include "PottsMeshGenerator.hpp"
#include "VolumeConstraintPottsUpdateRule.hpp"
#include "SurfaceAreaConstraintPottsUpdateRule.hpp"
#include "AdhesionPottsUpdateRule.hpp"
#include "DifferentialAdhesionPottsUpdateRule.hpp"

#include "CellIdWriter.hpp"
#include "CellMutationStatesWriter.hpp"

#include "Debug.hpp"

#include "CommandLineArguments.hpp"
#include "CommandLineArgumentsMocker.hpp"

#include "PetscSetupAndFinalize.hpp"

static const bool M_USING_COMMAND_LINE_ARGS = true;
static const double M_TIME_TO_STEADY_STATE = 10; //10
static const double M_TIME_FOR_SIMULATION = 100; //100
static const double M_NUM_CELLS_ACROSS = 20; //20 // this ^2 cells

class TestCellSorting: public AbstractCellBasedWithTimingsTestSuite
{
private:

    void RandomlyLabelCells(std::list<CellPtr>& rCells, boost::shared_ptr<AbstractCellProperty> pLabel, double labelledRatio)
    {
        for (std::list<CellPtr>::iterator cell_iter = rCells.begin();
             cell_iter != rCells.end();
             ++cell_iter)
        {
            if (RandomNumberGenerator::Instance()->ranf() < labelledRatio)
            {
               (*cell_iter)->AddCellProperty(pLabel);
            }
        }
    }

public:

   void NoTestCaBasedMonolayerCellSorting() throw (Exception)
    {
        double sim_index = 0;
        double cell_fluctuation = 1.0;
        if (M_USING_COMMAND_LINE_ARGS)
        {
          sim_index = (double) atof(CommandLineArguments::Instance()->GetStringCorrespondingToOption("-sim_index").c_str());
          cell_fluctuation = (double) atof(CommandLineArguments::Instance()->GetStringCorrespondingToOption("-noise").c_str());
        }
        RandomNumberGenerator::Instance()->Reseed(100.0*sim_index);

        //Create output directory
        std::stringstream out;
        out << sim_index << "_noise_" << cell_fluctuation;
        std::string output_directory = "CellSorting/Ca/" +  out.str();

        // Create a simple 2D PottsMesh
        unsigned domain_wide = 2*M_NUM_CELLS_ACROSS;

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
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<StochasticDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, location_indices.size(), p_differentiated_type);

        // Create cell population
        CaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);

        // Set population to output all data to results files
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();
        cell_population.AddPopulationWriter<HeterotypicBoundaryLengthWriter>();
        cell_population.AddPopulationWriter<CellPopulationAdjacencyMatrixWriter>();

        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory(output_directory);
        simulator.SetDt(0.01);
        simulator.SetSamplingTimestepMultiple(100);
        simulator.SetEndTime(M_TIME_TO_STEADY_STATE);

        // Add Division Rule
        boost::shared_ptr<AbstractCaBasedDivisionRule<2> > p_division_rule(new ShovingCaBasedDivisionRule<2>());
        cell_population.SetCaBasedDivisionRule(p_division_rule);

        // Add switching Update Rule
        MAKE_PTR(DifferentialAdhesionCaSwitchingUpdateRule<2u>, p_switching_update_rule);
        p_switching_update_rule->SetLabelledCellLabelledCellAdhesionEnergyParameter(0.1);
        p_switching_update_rule->SetLabelledCellCellAdhesionEnergyParameter(0.2);
        p_switching_update_rule->SetCellCellAdhesionEnergyParameter(0.1);
        p_switching_update_rule->SetCellBoundaryAdhesionEnergyParameter(0.2);
        p_switching_update_rule->SetLabelledCellBoundaryAdhesionEnergyParameter(0.4);
        p_switching_update_rule->SetTemperature(0.1);
        simulator.AddCaSwitchingUpdateRule(p_switching_update_rule);

        simulator.Solve();

        // Now label some cells
        boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<CellLabel>());
        RandomlyLabelCells(simulator.rGetCellPopulation().rGetCells(), p_state, 0.5);

        // modify parameters
        p_switching_update_rule->SetTemperature(0.1*cell_fluctuation);


        // Run simulation
        simulator.SetEndTime(M_TIME_TO_STEADY_STATE + M_TIME_FOR_SIMULATION);
        simulator.Solve();

        // Check that the same number of cells
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), M_NUM_CELLS_ACROSS*M_NUM_CELLS_ACROSS);

        // Test no births or deaths
        TS_ASSERT_EQUALS(simulator.GetNumBirths(), 0u);
        TS_ASSERT_EQUALS(simulator.GetNumDeaths(), 0u);
    }

    /**
     * Simulate a population of cells exhibiting cell sorting using the
     * Potts model.
     */
   void NoTestPottsMonolayerCellSorting() throw (Exception)
    {
        double sim_index = 0;
        double cell_fluctuation = 1.0;
        if (M_USING_COMMAND_LINE_ARGS)
        {
          sim_index = (double) atof(CommandLineArguments::Instance()->GetStringCorrespondingToOption("-sim_index").c_str());
          cell_fluctuation = (double) atof(CommandLineArguments::Instance()->GetStringCorrespondingToOption("-noise").c_str());
        }
        RandomNumberGenerator::Instance()->Reseed(100.0*sim_index);

        //Create output directory
        std::stringstream out;
        out << sim_index << "_noise_" << cell_fluctuation;
        std::string output_directory = "CellSorting/Potts/" +  out.str();

        // Create a simple 2D PottsMesh
        unsigned element_size = 4;
        unsigned domain_size = M_NUM_CELLS_ACROSS * element_size * 3; // Three times the initial domain size
        PottsMeshGenerator<2> generator(domain_size, M_NUM_CELLS_ACROSS, element_size, domain_size, M_NUM_CELLS_ACROSS, element_size);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<StochasticDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_differentiated_type);

        // Create cell population
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();
        cell_population.AddPopulationWriter<HeterotypicBoundaryLengthWriter>();
        cell_population.AddPopulationWriter<CellPopulationAdjacencyMatrixWriter>();

        // Set the Temperature
        cell_population.SetTemperature(0.2); //Default is 0.1

        // Set up cell-based simulation and output directory
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory(output_directory);

        // Set time step and end time for simulation
        simulator.SetDt(0.01);
        simulator.SetSamplingTimestepMultiple(100);
        simulator.SetEndTime(M_TIME_TO_STEADY_STATE);

        // Create update rules and pass to the simulation
        MAKE_PTR(VolumeConstraintPottsUpdateRule<2>, p_volume_constraint_update_rule);
        p_volume_constraint_update_rule->SetMatureCellTargetVolume(16); // i.e 4x4 cells
        p_volume_constraint_update_rule->SetDeformationEnergyParameter(0.1);
        simulator.AddPottsUpdateRule(p_volume_constraint_update_rule);

        MAKE_PTR(SurfaceAreaConstraintPottsUpdateRule<2>, p_surface_constraint_update_rule);
        p_surface_constraint_update_rule->SetMatureCellTargetSurfaceArea(16); // i.e 4x4 cells
        p_surface_constraint_update_rule->SetDeformationEnergyParameter(0.01);//0.01
        simulator.AddPottsUpdateRule(p_surface_constraint_update_rule);

        MAKE_PTR(DifferentialAdhesionPottsUpdateRule<2>, p_differential_adhesion_update_rule);
        p_differential_adhesion_update_rule->SetLabelledCellLabelledCellAdhesionEnergyParameter(0.1);
        p_differential_adhesion_update_rule->SetLabelledCellCellAdhesionEnergyParameter(0.5); // 1.0
        p_differential_adhesion_update_rule->SetCellCellAdhesionEnergyParameter(0.1); //0.1
        p_differential_adhesion_update_rule->SetCellBoundaryAdhesionEnergyParameter(0.2); // 1.0
        p_differential_adhesion_update_rule->SetLabelledCellBoundaryAdhesionEnergyParameter(1.0); // 2.0
        simulator.AddPottsUpdateRule(p_differential_adhesion_update_rule);

        // Run simulation
        simulator.Solve();

        // Now label some cells
        boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<CellLabel>());
        RandomlyLabelCells(simulator.rGetCellPopulation().rGetCells(), p_state, 0.5);

        // Adjust Parameters
        dynamic_cast <PottsBasedCellPopulation<2>*>(&(simulator.rGetCellPopulation()))->SetTemperature(0.2*cell_fluctuation);

        // Run simulation
        simulator.SetEndTime(M_TIME_TO_STEADY_STATE + M_TIME_FOR_SIMULATION);
        simulator.Solve();

        // Check that the same number of cells
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), M_NUM_CELLS_ACROSS*M_NUM_CELLS_ACROSS);

        // Test no births or deaths
        TS_ASSERT_EQUALS(simulator.GetNumBirths(), 0u);
        TS_ASSERT_EQUALS(simulator.GetNumDeaths(), 0u);
    }

   void NoTestNodeBasedMonolayerCellSorting() throw (Exception)
    {
        double sim_index = 0;
        double cell_fluctuation = 1.0;
        if (M_USING_COMMAND_LINE_ARGS)
        {
          sim_index = (double) atof(CommandLineArguments::Instance()->GetStringCorrespondingToOption("-sim_index").c_str());
          cell_fluctuation = (double) atof(CommandLineArguments::Instance()->GetStringCorrespondingToOption("-noise").c_str());
        }
        RandomNumberGenerator::Instance()->Reseed(100.0*sim_index);

        double cut_off_length = 2.5; //Extended to allow sorting for longer distances

        //Create output directory
        std::stringstream out;
        out << sim_index << "_noise_" << cell_fluctuation;
        std::string output_directory = "CellSorting/Node/" +  out.str();

        // Create a simple mesh
        HoneycombMeshGenerator generator(M_NUM_CELLS_ACROSS, M_NUM_CELLS_ACROSS, 0);
        TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh, cut_off_length);

        // Set up cells, one for each Node
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<StochasticDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);

        // Create cell population
        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        // Set population to output all data to results files
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();
        cell_population.AddPopulationWriter<HeterotypicBoundaryLengthWriter>();
        cell_population.AddPopulationWriter<CellPopulationAdjacencyMatrixWriter>();

        // Set up cell-based simulation and output directory
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory(output_directory);

        // Set time step and end time for simulation
        simulator.SetDt(1.0/200.0);
        simulator.SetSamplingTimestepMultiple(200);
        simulator.SetEndTime(M_TIME_TO_STEADY_STATE);

        // Create a force law and pass it to the simulation
        MAKE_PTR(DifferentialAdhesionGeneralisedLinearSpringForce<2>, p_differential_adhesion_force);
        p_differential_adhesion_force->SetMeinekeSpringStiffness(50.0);
        p_differential_adhesion_force->SetHomotypicLabelledSpringConstantMultiplier(1.0);
        p_differential_adhesion_force->SetHeterotypicSpringConstantMultiplier(0.1);
        p_differential_adhesion_force->SetCutOffLength(cut_off_length);
        simulator.AddForce(p_differential_adhesion_force);

        // Add some noise to avoid local minimum
        MAKE_PTR(RandomMotionForce<2>, p_random_force);
        p_random_force->SetMovementParameter(0.05); //0.1 causes dissasociation, 0.001 is not enough
        simulator.AddForce(p_random_force);

        // Run simulation
        simulator.Solve();

        // Now label some cells
        boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<CellLabel>());
        RandomlyLabelCells(simulator.rGetCellPopulation().rGetCells(), p_state, 0.5);

        // Adjust parameters
        p_random_force->SetMovementParameter(0.05*cell_fluctuation); //0.1 causes dissasociation


        // Run simulation
        simulator.SetEndTime(M_TIME_TO_STEADY_STATE + M_TIME_FOR_SIMULATION);
        simulator.Solve();

        // Check that the same number of cells
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), M_NUM_CELLS_ACROSS*M_NUM_CELLS_ACROSS);

        // Test no births or deaths
        TS_ASSERT_EQUALS(simulator.GetNumBirths(), 0u);
        TS_ASSERT_EQUALS(simulator.GetNumDeaths(), 0u);
    }


//   void NoTestMeshBasedMolayerCellSorting() throw (Exception)
//    {
//        double sim_index = 0;
//        double cell_fluctuation = 1.0;
//        if (M_USING_COMMAND_LINE_ARGS)
//        {
//          sim_index = (double) atof(CommandLineArguments::Instance()->GetStringCorrespondingToOption("-sim_index").c_str());
//          cell_fluctuation = (double) atof(CommandLineArguments::Instance()->GetStringCorrespondingToOption("-noise").c_str());
//        }
//        RandomNumberGenerator::Instance()->Reseed(100.0*sim_index);
//
//        //Create output directory
//        std::stringstream out;
//        out << sim_index << "_noise_" << cell_fluctuation;
//        std::string output_directory = "CellSorting/Mesh/" +  out.str();
//
//        // Create a simple mesh
//        HoneycombMeshGenerator generator(M_NUM_CELLS_ACROSS, M_NUM_CELLS_ACROSS, 0);
//        MutableMesh<2,2>* p_mesh = generator.GetMesh();
//
//        std::vector<CellPtr> cells;
//        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
//        CellsGenerator<StochasticDurationGenerationBasedCellCycleModel, 2> cells_generator;
//        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), p_differentiated_type);
//
//        // Create cell population
//        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);
//        cell_population.AddPopulationWriter<CellPopulationAdjacencyMatrixWriter>();
//
//        // Set population to output all data to results files
//        cell_population.AddCellWriter<CellIdWriter>();
//        cell_population.AddCellWriter<CellMutationStatesWriter>();
//        cell_population.AddPopulationWriter<HeterotypicBoundaryLengthWriter>();
//
//        // Set up cell-based simulation and output directory
//        OffLatticeSimulation<2> simulator(cell_population);
//        simulator.SetOutputDirectory(output_directory);
//
//        // Set time step and end time for simulation
//        simulator.SetDt(1.0/200.0);
//        simulator.SetSamplingTimestepMultiple(200);
//        simulator.SetEndTime(M_TIME_TO_STEADY_STATE);
//
//        // Create a force law and pass it to the simulation
//        MAKE_PTR(DifferentialAdhesionGeneralisedLinearSpringForce<2>, p_differential_adhesion_force);
//        p_differential_adhesion_force->SetMeinekeSpringStiffness(50.0);
//        p_differential_adhesion_force->SetHomotypicLabelledSpringConstantMultiplier(1.0);
//        p_differential_adhesion_force->SetHeterotypicSpringConstantMultiplier(0.1);
//        simulator.AddForce(p_differential_adhesion_force);
//
//        // Add some noise to avoid local minimum
//        MAKE_PTR(RandomMotionForce<2>, p_random_force);
//        p_random_force->SetMovementParameter(0.02);
//        simulator.AddForce(p_random_force);
//
//        // Run simulation
//        simulator.Solve();
//
//        // Now label some cells
//        boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<CellLabel>());
//        RandomlyLabelCells(simulator.rGetCellPopulation().rGetCells(), p_state, 0.5);
//
//        // Adjust parameters
//        p_random_force->SetMovementParameter(0.02*cell_fluctuation); //0.1 causes dissasociation
//
//        // Run simulation
//        simulator.SetEndTime(M_TIME_TO_STEADY_STATE + M_TIME_FOR_SIMULATION);
//        simulator.Solve();
//
//        // Check that the same number of cells
//        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), M_NUM_CELLS_ACROSS*M_NUM_CELLS_ACROSS);
//
//        // Test no births or deaths
//        TS_ASSERT_EQUALS(simulator.GetNumBirths(), 0u);
//        TS_ASSERT_EQUALS(simulator.GetNumDeaths(), 0u);
//    }

   void noTestMeshBasedWithGhostsMonolayerCellSorting() throw (Exception)
    {
        double sim_index = 0;
        double cell_fluctuation = 1.0;
        if (M_USING_COMMAND_LINE_ARGS)
        {
          sim_index = (double) atof(CommandLineArguments::Instance()->GetStringCorrespondingToOption("-sim_index").c_str());
          cell_fluctuation = (double) atof(CommandLineArguments::Instance()->GetStringCorrespondingToOption("-noise").c_str());
        }
        RandomNumberGenerator::Instance()->Reseed(100.0*sim_index);

        //Create output directory
        std::stringstream out;
        out << sim_index << "_noise_" << cell_fluctuation;
        std::string output_directory = "CellSorting/MeshGhost/" +  out.str();

        // Create a simple mesh
        unsigned num_ghosts = 20;
        HoneycombMeshGenerator generator(M_NUM_CELLS_ACROSS, M_NUM_CELLS_ACROSS, num_ghosts);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        // Set up cells, one for each non ghost Node
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<StochasticDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, location_indices.size(), p_differentiated_type);

        // Create cell population
        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);
        cell_population.AddPopulationWriter<CellPopulationAdjacencyMatrixWriter>();

        // Set population to output all data to results files
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();
        cell_population.AddPopulationWriter<HeterotypicBoundaryLengthWriter>();

        // Set up cell-based simulation and output directory
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory(output_directory);

        // Set time step and end time for simulation
        simulator.SetDt(1.0/200.0);
        simulator.SetSamplingTimestepMultiple(200);
        simulator.SetEndTime(M_TIME_TO_STEADY_STATE);

        // Create a force law and pass it to the simulation
        MAKE_PTR(DifferentialAdhesionGeneralisedLinearSpringForce<2>, p_differential_adhesion_force);
        p_differential_adhesion_force->SetMeinekeSpringStiffness(50.0);
        p_differential_adhesion_force->SetHomotypicLabelledSpringConstantMultiplier(1.0);
        p_differential_adhesion_force->SetHeterotypicSpringConstantMultiplier(0.1);
        simulator.AddForce(p_differential_adhesion_force);

        // Add some noise to avoid local minimum
        MAKE_PTR(RandomMotionForce<2>, p_random_force);
        p_random_force->SetMovementParameter(0.05);
        simulator.AddForce(p_random_force);

        // Run simulation
        simulator.Solve();

        // Now label some cells
        boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<CellLabel>());
        RandomlyLabelCells(simulator.rGetCellPopulation().rGetCells(), p_state, 0.5);

        // Adjust parameters
        p_random_force->SetMovementParameter(0.05*cell_fluctuation); //0.1 causes dissasociation

        // Run simulation
        simulator.SetEndTime(M_TIME_TO_STEADY_STATE + M_TIME_FOR_SIMULATION);
        simulator.Solve();

        // Check that the same number of cells
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), M_NUM_CELLS_ACROSS*M_NUM_CELLS_ACROSS);

        // Test no births or deaths
        TS_ASSERT_EQUALS(simulator.GetNumBirths(), 0u);
        TS_ASSERT_EQUALS(simulator.GetNumDeaths(), 0u);
    }


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
        double sim_index = 0;
        double cell_fluctuation = 1.0;
        if (M_USING_COMMAND_LINE_ARGS)
        {
          sim_index = (double) atof(CommandLineArguments::Instance()->GetStringCorrespondingToOption("-sim_index").c_str());
          cell_fluctuation = (double) atof(CommandLineArguments::Instance()->GetStringCorrespondingToOption("-noise").c_str());
        }
        RandomNumberGenerator::Instance()->Reseed(100.0*sim_index);

        //Create output directory
        std::stringstream out;
        out << sim_index << "_noise_" << cell_fluctuation;
        std::string output_directory = "CellSorting/Vertex/" +  out.str();

        // Create a simple 2D MutableVertexMesh
        HoneycombVertexMeshGenerator generator(M_NUM_CELLS_ACROSS, M_NUM_CELLS_ACROSS);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();
        p_mesh->SetCellRearrangementThreshold(0.1);

        // Slows things down but can use a larger timestep and diffusion forces
        //p_mesh->SetCheckForInternalIntersections(true);

        // Set up cells, one for each VertexElement
        std::vector<CellPtr> cells;
        boost::shared_ptr<AbstractCellProperty> p_cell_type(CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>());
        CellsGenerator<StochasticDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_cell_type);

        for (unsigned i=0; i<cells.size(); i++)
        {
            // Set a target area rather than setting a growth modifier. (the modifiers don't work correctly as making very long G1 phases)
            cells[i]->GetCellData()->SetItem("target area", 1.0);
        }

        // Create cell population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set population to output all data to results files
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();
        cell_population.AddPopulationWriter<HeterotypicBoundaryLengthWriter>();
        cell_population.AddPopulationWriter<CellPopulationAdjacencyMatrixWriter>();

        // Set up cell-based simulation and output directory
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory(output_directory);

        // Set time step and end time for simulation
        simulator.SetDt(1.0/200.0);
        simulator.SetSamplingTimestepMultiple(200);
        simulator.SetEndTime(M_TIME_TO_STEADY_STATE);

        // Set up force law and pass it to the simulation
        MAKE_PTR(NagaiHondaDifferentialAdhesionForce<2>, p_force);
        p_force->SetNagaiHondaDeformationEnergyParameter(50.0);
        p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(1.0);
        p_force->SetNagaiHondaCellCellAdhesionEnergyParameter(1.0);
        p_force->SetNagaiHondaLabelledCellCellAdhesionEnergyParameter(2.0);
        p_force->SetNagaiHondaLabelledCellLabelledCellAdhesionEnergyParameter(1.0);
        p_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(10.0);
        p_force->SetNagaiHondaLabelledCellBoundaryAdhesionEnergyParameter(20.0);
        simulator.AddForce(p_force);

        // Add some noise to avoid local minimum
        MAKE_PTR(RandomMotionForce<2>, p_random_force);
        p_random_force->SetMovementParameter(0.05);
        simulator.AddForce(p_random_force);

        // Run simulation
        simulator.Solve();

        // Now label some cells
        boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<CellLabel>());
        RandomlyLabelCells(simulator.rGetCellPopulation().rGetCells(), p_state, 0.5);

        // Adjust parameters
        p_random_force->SetMovementParameter(0.05*cell_fluctuation);

        // Run simulation
        simulator.SetEndTime(M_TIME_TO_STEADY_STATE + M_TIME_FOR_SIMULATION);
        simulator.Solve();


        // Check that the same number of cells
        TS_ASSERT_EQUALS(simulator.rGetCellPopulation().GetNumRealCells(), M_NUM_CELLS_ACROSS*M_NUM_CELLS_ACROSS);

        // Test no births or deaths
        TS_ASSERT_EQUALS(simulator.GetNumBirths(), 0u);
        TS_ASSERT_EQUALS(simulator.GetNumDeaths(), 0u);
   }
};
