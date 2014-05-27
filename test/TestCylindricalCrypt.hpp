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
#include "ContactInhibitionGenerationBasedCellCycleModel.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "NagaiHondaForce.hpp"
#include "ShortAxisStemHorizontalVertexBasedDivisionRule.hpp"

#include "Cylindrical2dNodesOnlyMesh.hpp"
#include "NodesOnlyMesh.hpp"

#include "MeshBasedCellPopulationWithGhostNodes.hpp"

#include "PeriodicNodeBasedBoundaryCondition.hpp"
#include "StemCellRetainerForce.hpp"

#include "PottsMeshGenerator.hpp"
#include "SimpleWntCellCycleModel.hpp"
#include "OnLatticeSimulation.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "VolumeConstraintPottsUpdateRule.hpp"
#include "AdhesionPottsUpdateRule.hpp"

#include "SimpleTargetAreaModifier.hpp"
#include "VolumeTrackingModifier.hpp"
#include "DiffusionCaUpdateRule.hpp"

#include "CellVolumesWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "CellIdWriter.hpp"


#include "PlaneBasedCellKiller.hpp"
#include "PlaneBoundaryCondition.hpp"
#include "StemCellRetainerPottsUpdateRule.hpp"

#include "PetscSetupAndFinalize.hpp"

#include "Warnings.hpp"
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

    void GenerateStemCells(unsigned num_cells, std::vector<CellPtr>& rCells, double EquilibriumVolume)
	{
    	double typical_stem_cell_cycle_duration = 24.0;

    	boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());
    	boost::shared_ptr<AbstractCellProperty> p_cell_type(CellPropertyRegistry::Instance()->Get<StemCellProliferativeType>());

		for (unsigned i=0; i<num_cells; i++)
		{
			ContactInhibitionGenerationBasedCellCycleModel* p_model = new ContactInhibitionGenerationBasedCellCycleModel();
			p_model->SetDimension(2);
			p_model->SetMaxTransitGenerations(3);

			p_model->SetEquilibriumVolume(EquilibriumVolume);
			// 0.1 -> No CI!!!!
			p_model->SetQuiescentVolumeFraction(0.1); //0.8

			CellPtr p_cell(new Cell(p_state, p_model));
			p_cell->SetCellProliferativeType(p_cell_type);
			double birth_time = - RandomNumberGenerator::Instance()->ranf() * typical_stem_cell_cycle_duration;
			p_cell->SetBirthTime(birth_time);
			rCells.push_back(p_cell);
		}
	}

public:

    const double mEndTime = 1200.0;

    void TestVertexCrypt() throw (Exception)
    {
        double crypt_length = 12;
        double crypt_width = 5.0;

        // Create mesh
        unsigned cells_across = 6;

        CylindricalHoneycombVertexMeshGenerator generator(cells_across, 1, true);
        Cylindrical2dVertexMesh* p_mesh = generator.GetCylindricalMesh();

        p_mesh->Scale(crypt_width/cells_across,1.0);

        // Create cells
        std::vector<CellPtr> cells;
        GenerateStemCells(cells_across,cells,0.9); //1.0

        // Create tissue
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddPopulationWriter<CellProliferativeTypesCountWriter>();
		cell_population.AddCellWriter<CellVolumesWriter>();
		cell_population.AddCellWriter<CellIdWriter>();

        // Set the division rule so stem cells divide horizontally
        MAKE_PTR(ShortAxisStemHorizontalVertexBasedDivisionRule<2>,p_division_rule);
        cell_population.SetVertexBasedDivisionRule(p_division_rule);

        // Create crypt simulation from cell population
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetDt(1.0/500.0);
        simulator.SetSamplingTimestepMultiple(500);
        simulator.SetEndTime(mEndTime);
        simulator.SetOutputDirectory("CylindricalCrypt/Vertex");
        simulator.SetOutputDivisionLocations(true);
        simulator.SetOutputCellVelocities(true);

        //Add Volume Tracking Moddifer
        MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        // Create a force law and pass it to the simulation
        MAKE_PTR(NagaiHondaForce<2>, p_nagai_honda_force);
        simulator.AddForce(p_nagai_honda_force);

        // A NagaiHondaForce has to be used together with an AbstractTargetAreaModifier
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        // Create a force law to retain stem cells niche and pass it to the simulation
        MAKE_PTR(StemCellRetainerForce<2>, p_retainer_force);
        p_retainer_force->SetStemCellForceMagnitudeParameter(1.0);
        simulator.AddForce(p_retainer_force);

        // Solid Base Boundary Condition
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

    void TestMeshBasedCrypt() throw (Exception)
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
        GenerateStemCells(cells_across,cells,sqrt(3.0)/2.0);

        // Create tissue
        MeshBasedCellPopulationWithGhostNodes<2> cell_population(*p_mesh, cells, location_indices);
        cell_population.AddPopulationWriter<CellProliferativeTypesCountWriter>();
		cell_population.AddCellWriter<CellVolumesWriter>();
		cell_population.AddCellWriter<CellIdWriter>();

        // Create simulation from cell population
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetDt(1.0/200.0);
        simulator.SetSamplingTimestepMultiple(200);
        simulator.SetEndTime(mEndTime);
        simulator.SetOutputDirectory("CylindricalCrypt/Mesh");
        simulator.SetOutputDivisionLocations(true);
        simulator.SetOutputCellVelocities(true);

        //Add Volume Tracking Moddifer
        MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        // Create a force law and pass it to the simulation
        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetMeinekeSpringStiffness(50.0);
        simulator.AddForce(p_linear_force);

        // Create a force law to retain stem cells niche and pass it to the simulation
        MAKE_PTR(StemCellRetainerForce<2>, p_retainer_force);
        p_retainer_force->SetStemCellForceMagnitudeParameter(50.0);
        simulator.AddForce(p_retainer_force);

        // Solid Base Boundary Condition
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bcs, (&cell_population, zero_vector<double>(2), -unit_vector<double>(2,1)));
        p_bcs->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(p_bcs);

        // Sloughing killer
        MAKE_PTR_ARGS(PlaneBasedCellKiller<2>, p_killer, (&cell_population, crypt_length*unit_vector<double>(2,1), unit_vector<double>(2,1)));
        simulator.AddCellKiller(p_killer);

        // Run simulation
        simulator.Solve();
    }


    void TestNodeBasedCrypt() throw (Exception)
    {
        double crypt_length = 12-0.5;
        double crypt_width = 5.0;

        // Create a simple mesh
        unsigned cells_across = 6;
        HoneycombMeshGenerator generator(cells_across, 1, 0, crypt_width/cells_across);
        TetrahedralMesh<2,2>* p_generating_mesh = generator.GetMesh();

        // Convert this to a Cylindrical2dNodesOnlyMesh
        Cylindrical2dNodesOnlyMesh* p_mesh = new Cylindrical2dNodesOnlyMesh(crypt_width);
        p_mesh->ConstructNodesWithoutMesh(*p_generating_mesh,crypt_width/2.);

        // Get location indices corresponding to real cells
        std::vector<unsigned> location_indices = generator.GetCellLocationIndices();

        // Create cells
        std::vector<CellPtr> cells;
        GenerateStemCells(cells_across,cells,0.9*M_PI*0.25); // r=0.5 remove 0.9*!!!!!

        // Create a node-based cell population
        NodeBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddPopulationWriter<CellProliferativeTypesCountWriter>();
		cell_population.AddCellWriter<CellVolumesWriter>();
		cell_population.AddCellWriter<CellIdWriter>();

        for (unsigned index = 0; index < cell_population.rGetMesh().GetNumNodes(); index++)
        {
            //PRINT_VARIABLE(index);
        	cell_population.rGetMesh().GetNode(index)->SetRadius(0.5);
        }

        // Create simulation from cell population
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetDt(1.0/200.0);
        simulator.SetSamplingTimestepMultiple(200);
        simulator.SetEndTime(mEndTime);
        simulator.SetOutputDirectory("CylindricalCrypt/Node");
        simulator.SetOutputDivisionLocations(true);
        simulator.SetOutputCellVelocities(true);

        //Add Volume Tracking Moddifer
        MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        // Create a force law and pass it to the simulation
        MAKE_PTR(RepulsionForce<2>, p_repulsion_force);
        p_repulsion_force->SetMeinekeSpringStiffness(50.0);
        simulator.AddForce(p_repulsion_force);

        // Create a force law to retain stem cells niche and pass it to the simulation
        MAKE_PTR(StemCellRetainerForce<2>, p_retainer_force);
        p_retainer_force->SetStemCellForceMagnitudeParameter(50.0);
        simulator.AddForce(p_retainer_force);

        // Solid Base Boundary Condition
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bcs, (&cell_population, zero_vector<double>(2), -unit_vector<double>(2,1)));
        p_bcs->SetUseJiggledNodesOnPlane(true);
        simulator.AddCellPopulationBoundaryCondition(p_bcs);

        // Sloughing killer
        MAKE_PTR_ARGS(PlaneBasedCellKiller<2>, p_killer, (&cell_population, crypt_length*unit_vector<double>(2,1), unit_vector<double>(2,1)));
        simulator.AddCellKiller(p_killer);

        // Periodic BCS to move nodes round periodic boundary
        MAKE_PTR_ARGS(PeriodicNodeBasedBoundaryCondition, p_periodic_bcs,(&cell_population, crypt_width));
        simulator.AddCellPopulationBoundaryCondition(p_periodic_bcs);

        // Run simulation
        simulator.Solve();

        // Clear memory
        delete p_mesh;
    }

    void TestPottsCrypt() throw (Exception)
    {
        double crypt_length = 12;
        double crypt_width = 6;
        unsigned cells_across = 6;

        // Create a simple 2D PottsMesh (periodic in x)
        PottsMeshGenerator<2> generator(crypt_width*4, cells_across, 4, (crypt_length+2)*4, 1, 4, 1, 1, 1, true);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        GenerateStemCells(cells_across,cells,18.0); //16.0

        // Create cell population
        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetNumSweepsPerTimestep(10);
        cell_population.SetTemperature(0.001);
        cell_population.AddPopulationWriter<CellProliferativeTypesCountWriter>();
		cell_population.AddCellWriter<CellVolumesWriter>();
		cell_population.AddCellWriter<CellIdWriter>();

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("CylindricalCrypt/Potts");
        simulator.SetDt(0.1);
        simulator.SetSamplingTimestepMultiple(10);
        simulator.SetEndTime(mEndTime);
        simulator.SetOutputDivisionLocations(true);
        simulator.SetOutputCellVelocities(true);

        //Add Volume Tracking Moddifer
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

        // Restrain Stem Cells at Base of Crypt
        MAKE_PTR(StemCellRetainerPottsUpdateRule<2>, p_stem_cell_retainer_update_rule);
        p_stem_cell_retainer_update_rule->SetStemCellRestraintParameter(1000.0);
        simulator.AddPottsUpdateRule(p_stem_cell_retainer_update_rule);

        // Run simulation
        simulator.Solve();
    }



    void TestCaCrypt() throw (Exception)
    {
        unsigned cells_across = 6;
        unsigned crypt_length = 12;

        // Create a simple 2D PottsMesh (periodic in x)
        PottsMeshGenerator<2> generator(6, 0, 0, crypt_length+1, 0, 0, 1, 0, 0, true);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        GenerateStemCells(cells_across,cells,1.0);

        // Specify where cells lie
        std::vector<unsigned> location_indices;
        for (unsigned index=0; index<cells_across; index++)
        {
            location_indices.push_back(index);
        }

        // Create cell population
        CaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);
        cell_population.AddPopulationWriter<CellProliferativeTypesCountWriter>();
		cell_population.AddCellWriter<CellVolumesWriter>();
		cell_population.AddCellWriter<CellIdWriter>();

        // Set up cell-based simulation
        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("CylindricalCrypt/Ca");
        simulator.SetDt(0.1);
        simulator.SetSamplingTimestepMultiple(10);
        simulator.SetEndTime(mEndTime);
        simulator.SetOutputDivisionLocations(true);
        simulator.SetOutputCellVelocities(true);

        //Add Volume Tracking Moddifer
        MAKE_PTR(VolumeTrackingModifier<2>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        // Add update rule
        MAKE_PTR(DiffusionCaUpdateRule<2>, p_diffusion_update_rule);
        p_diffusion_update_rule->SetDiffusionParameter(0.1);
        simulator.AddCaUpdateRule(p_diffusion_update_rule);

        // Sloughing killer
        MAKE_PTR_ARGS(PlaneBasedCellKiller<2>, p_killer, (&cell_population, crypt_length*unit_vector<double>(2,1), unit_vector<double>(2,1)));
        simulator.AddCellKiller(p_killer);

        // Run simulation
        simulator.Solve();
    }
};

#endif /*TESTCRYPTSIMULATIONCOMPARISON_HPP_*/
