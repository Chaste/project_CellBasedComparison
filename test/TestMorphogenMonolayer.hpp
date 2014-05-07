#include <cxxtest/TestSuite.h>
#include "CellBasedSimulationArchiver.hpp"
#include "SmartPointers.hpp"
#include "AbstractCellBasedTestSuite.hpp"

#include "TransitCellProliferativeType.hpp"

#include "CellIdWriter.hpp"
#include "VoronoiDataWriter.hpp"
#include "CellMutationStatesWriter.hpp"

#include "NeumannPDEModifier.hpp"
#include "CellwiseSourceMorphogenPde.hpp"
#include "VolumeTrackingModifier.hpp"

#include "MorphogenDependentCellCycleModel.hpp"
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

#include "MultipleCaBasedCellPopulation.hpp"
#include "DiffusionMultipleCaUpdateRule.hpp"
#include "CellMutationStatesCountWriter.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestMorphogenMonolayers : public AbstractCellBasedTestSuite
{

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
    	MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        boost::shared_ptr<AbstractCellProperty> p_label(CellPropertyRegistry::Instance()->Get<CellLabel>());

        for (unsigned i=0; i<num_cells; i++)
        {
        	double typical_cell_cycle_duration = 12.0;

        	MorphogenDependentCellCycleModel* p_cycle_model = new MorphogenDependentCellCycleModel();
            p_cycle_model->SetDimension(2);
            p_cycle_model->SetThresholdMorphogen(0.0);
            p_cycle_model->SetEquilibriumVolume(16.0);
            p_cycle_model->SetQuiescentVolumeFraction(0.5);

            CellPtr p_cell(new Cell(p_state, p_cycle_model));
            p_cell->SetCellProliferativeType(p_transit_type);

            double birth_time = - RandomNumberGenerator::Instance()->ranf() * typical_cell_cycle_duration;
 			p_cell->SetBirthTime(birth_time);

 			p_cell->InitialiseCellCycleModel();

 			if (i<(double)num_cells/2.0)
 			{
 				p_cell->AddCellProperty(p_label);
 			}

 			rCells.push_back(p_cell);
        }
 	}


	/*
	 * Simulates diffusion on a growing monolayer.
	 */

public:
    void noTestVertexBasedMonolayer() throw (Exception)
    {
        // Create Mesh
        HoneycombVertexMeshGenerator generator(10, 10);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        // Create Cells
        std::vector<CellPtr> cells;
        GenerateCells(p_mesh->GetNumElements(),cells);

        // Create Population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();

        // Create Simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("Monolayers/Vertex");
        simulator.SetDt(1.0/200.0);
        simulator.SetSamplingTimestepMultiple(200);
        simulator.SetEndTime(50.0);//20

        // Create Forces and pass to simulation NOTE: PARAMETERS CHOSEN TO GET CIRCULAR MONOLAYER
        MAKE_PTR(NagaiHondaForce<2>, p_force);
        p_force->SetNagaiHondaDeformationEnergyParameter(55.0);
		p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(0.0);
		p_force->SetNagaiHondaCellCellAdhesionEnergyParameter(6.0);
		p_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(12.0);
        simulator.AddForce(p_force);

        // Create Modifiers and pass to simulation

        // Create a pde modifier and pass it to the simulation Add this first so in place for SimpleTargetArea one (calls cell pop update)

        // Set up PDE and pass to simulation via modifier (uniform secretion at each labeled cell)
        CellwiseSourceMorphogenPde<2> pde(cell_population, 0.1);
        MAKE_PTR_ARGS(NeumannPDEModifier<2>, p_pde_modifier, (&pde));
        simulator.AddSimulationModifier(p_pde_modifier);

        // Add volume tracking modifier for CI cell cycle model.
        MAKE_PTR(VolumeTrackingModifier<2>, p_volume_tracking_modifier);
        simulator.AddSimulationModifier(p_volume_tracking_modifier);

        // A NagaiHondaForce has to be used together with an AbstractTargetAreaModifier
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);



        simulator.Solve();
    }

    void noTestNodeBasedMonolayer() throw (Exception)
    {
        HoneycombMeshGenerator generator(10, 10);
        MutableMesh<2,2>* p_generating_mesh = generator.GetMesh();
        NodesOnlyMesh<2>* p_mesh = new NodesOnlyMesh<2>;
        p_mesh->ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);

        std::vector<CellPtr> cells;
        GenerateCells(p_mesh->GetNumNodes(),cells);

        NodeBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddCellWriter<CellIdWriter>();
        //TODO should be CellMutationStatesWriter
        cell_population.AddPopulationWriter<CellMutationStatesCountWriter>();

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("Monolayers/Node");
        simulator.SetDt(1.0/120.0);
        simulator.SetSamplingTimestepMultiple(120);
        simulator.SetEndTime(50.0);

        MAKE_PTR(RepulsionForce<2>, p_force);
        simulator.AddForce(p_force);

        // Set up PDE and pass to simulation via modifier (uniform secretion at each labeled cell)
		CellwiseSourceMorphogenPde<2> pde(cell_population, 0.1);

		// Create a pde modifier and pass it to the simulation
		MAKE_PTR_ARGS(NeumannPDEModifier<2>, p_pde_modifier, (&pde));
		simulator.AddSimulationModifier(p_pde_modifier);

		// Add volume tracking modifier for CI cell cycle model.
		MAKE_PTR(VolumeTrackingModifier<2>, p_volume_tracking_modifier);
		simulator.AddSimulationModifier(p_volume_tracking_modifier);

        simulator.Solve();

        delete p_mesh; // to stop memory leaks
    }

    void noTestMeshBasedMonolayer() throw (Exception)
    {
        HoneycombMeshGenerator generator(10,10);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        GenerateCells(p_mesh->GetNumNodes(),cells);

        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

		// Set population to output all data to results files
		cell_population.AddCellWriter<CellIdWriter>();
		cell_population.AddCellWriter<CellMutationStatesWriter>();

        cell_population.SetWriteVtkAsPoints(false);
        cell_population.AddPopulationWriter<VoronoiDataWriter>();


        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("Monolayers/Mesh");
        simulator.SetDt(1.0/120.0);
        simulator.SetSamplingTimestepMultiple(120);
        simulator.SetEndTime(50.0);

        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        p_force->SetCutOffLength(1.5);
        simulator.AddForce(p_force);

        // Set up PDE and pass to simulation via modifier (uniform secretion at each labeled cell)
        CellwiseSourceMorphogenPde<2> pde(cell_population, 0.1);

        // Create a pde modifier and pass it to the simulation
        MAKE_PTR_ARGS(NeumannPDEModifier<2>, p_pde_modifier, (&pde));
        simulator.AddSimulationModifier(p_pde_modifier);

        // Add volume tracking modifier for CI cell cycle model.
        MAKE_PTR(VolumeTrackingModifier<2>, p_volume_tracking_modifier);
        simulator.AddSimulationModifier(p_volume_tracking_modifier);

        simulator.Solve();
    }

    void TestPottsBasedMonolayer() throw (Exception)
    {
        PottsMeshGenerator<2> generator(200, 10, 4, 200, 10, 4);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        GenerateCells(p_mesh->GetNumElements(),cells);

        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetTemperature(0.1);
        cell_population.SetNumSweepsPerTimestep(100);

		// Set population to output all data to results files
		cell_population.AddCellWriter<CellIdWriter>();
		cell_population.AddCellWriter<CellMutationStatesWriter>();

        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("Monolayers/Potts");
        simulator.SetDt(1.0);
        simulator.SetSamplingTimestepMultiple(1);
        simulator.SetEndTime(50.0);

        MAKE_PTR(VolumeConstraintPottsUpdateRule<2>, p_volume_constraint_update_rule);
        simulator.AddPottsUpdateRule(p_volume_constraint_update_rule);
        MAKE_PTR(SurfaceAreaConstraintPottsUpdateRule<2>, p_surface_area_update_rule);
        simulator.AddPottsUpdateRule(p_surface_area_update_rule);
        MAKE_PTR(AdhesionPottsUpdateRule<2>, p_adhesion_update_rule);
        simulator.AddPottsUpdateRule(p_adhesion_update_rule);

        // Set up PDE and pass to simulation via modifier (uniform secretion at each labeled cell)
        CellwiseSourceMorphogenPde<2> pde(cell_population, 0.1);

        // Create a pde modifier and pass it to the simulation
        MAKE_PTR_ARGS(NeumannPDEModifier<2>, p_pde_modifier, (&pde));
        simulator.AddSimulationModifier(p_pde_modifier);

        // Add volume tracking modifier for CI cell cycle model.
        MAKE_PTR(VolumeTrackingModifier<2>, p_volume_tracking_modifier);
        simulator.AddSimulationModifier(p_volume_tracking_modifier);


        simulator.Solve();
    }

    void noTestCaBasedMonolayer() throw (Exception)
    {
        // Create a simple 2D PottsMesh
    	unsigned domain_wide = 50;

    	PottsMeshGenerator<2> generator(domain_wide, 0, 0, domain_wide, 0, 0);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Specify where cells lie
        unsigned cells_wide = 10;
        std::vector<unsigned> location_indices;
        for (unsigned i=0; i<cells_wide; i++)
        {
        	for (unsigned j=0; j<cells_wide; j++)
			{
        		unsigned offset = (domain_wide+1) * (domain_wide-cells_wide)/2;
				location_indices.push_back(offset + j + i * domain_wide );
			}
        }

        MAKE_PTR(TransitCellProliferativeType, p_transit_type);

        std::vector<CellPtr> cells;
        GenerateCells(location_indices.size(),cells);

        // Create cell population
        MultipleCaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);

		// Set population to output all data to results files
		cell_population.AddCellWriter<CellIdWriter>();
		cell_population.AddCellWriter<CellMutationStatesWriter>();

        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("Monolayers/CA");
        simulator.SetDt(0.1);
        simulator.SetSamplingTimestepMultiple(10);
        simulator.SetEndTime(50.0);

        // Adding update rule(s).
        MAKE_PTR(DiffusionMultipleCaUpdateRule<2u>, p_diffusion_update_rule);
        p_diffusion_update_rule->SetDiffusionParameter(0.1);

        simulator.AddMultipleCaUpdateRule(p_diffusion_update_rule);

        // Set up PDE and pass to simulation via modifier (uniform secretion at each labeled cell)
	   CellwiseSourceMorphogenPde<2> pde(cell_population, 0.1);

	   // Create a pde modifier and pass it to the simulation
	   MAKE_PTR_ARGS(NeumannPDEModifier<2>, p_pde_modifier, (&pde));
	   simulator.AddSimulationModifier(p_pde_modifier);

	   // Add volume tracking modifier for CI cell cycle model.
	   MAKE_PTR(VolumeTrackingModifier<2>, p_volume_tracking_modifier);
	   simulator.AddSimulationModifier(p_volume_tracking_modifier);

        // Run simulation
        simulator.Solve();
    }

};
