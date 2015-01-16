#include <cxxtest/TestSuite.h>
#include "CellBasedSimulationArchiver.hpp"
#include "SmartPointers.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"

#include "TransitCellProliferativeType.hpp"

#include "CellIdWriter.hpp"
#include "VoronoiDataWriter.hpp"
#include "CellMutationStatesWriter.hpp"

#include "ParabolicPDEModifier.hpp"
#include "CellwiseSourceParabolicPde.hpp"
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

#include "CaBasedCellPopulation.hpp"
#include "DiffusionCaUpdateRule.hpp"

#include "PetscSetupAndFinalize.hpp"

static const double M_TIME_FOR_SIMULATION = 10.0;
static const double M_NUM_CELLS_ACROSS = 5; // this ^2 cells
static const double M_UPTAKE_RATE = 0.1;
static const double M_DIFFUSION_CONSTANT = 1.0;
static const double M_DUDT_COEFFICIENT = 10.0;

class TestParabolicMorphogenMonolayers : public AbstractCellBasedWithTimingsTestSuite
{
private:
    void GenerateCells(unsigned num_cells, std::vector<CellPtr>& rCells, double EquilibriumVolume)
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
            p_cycle_model->SetEquilibriumVolume(EquilibriumVolume);
            p_cycle_model->SetQuiescentVolumeFraction(0.5);

            CellPtr p_cell(new Cell(p_state, p_cycle_model));
            p_cell->SetCellProliferativeType(p_transit_type);

            double birth_time = - RandomNumberGenerator::Instance()->ranf() * typical_cell_cycle_duration;
             p_cell->SetBirthTime(birth_time);

             p_cell->InitialiseCellCycleModel();

             if (i<(double)num_cells/2.0)
             {
                 p_cell->AddCellProperty(p_label);
                 p_cell->GetCellData()->SetItem("morphogen",1.0);
             }
             else
             {
                 p_cell->GetCellData()->SetItem("morphogen",0.0);
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
        HoneycombVertexMeshGenerator generator(M_NUM_CELLS_ACROSS, M_NUM_CELLS_ACROSS);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        // Create Cells
        std::vector<CellPtr> cells;
        GenerateCells(p_mesh->GetNumElements(),cells, 1.0);

        // Create Population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();

        // Create Simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("ParabolicMonolayers/Vertex");
        simulator.SetDt(1.0/200.0);
        simulator.SetSamplingTimestepMultiple(200);
        simulator.SetEndTime(M_TIME_FOR_SIMULATION);//20

        // Create Forces and pass to simulation NOTE: PARAMETERS CHOSEN TO GET CIRCULAR MONOLAYER
        MAKE_PTR(NagaiHondaForce<2>, p_force);
        p_force->SetNagaiHondaDeformationEnergyParameter(55.0);
        p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(0.0);
        p_force->SetNagaiHondaCellCellAdhesionEnergyParameter(6.0);
        p_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(12.0);
        simulator.AddForce(p_force);

        // Create Modifiers and pass to simulation

        // Set up PDE and pass to simulation via modifier (uniform secretion at each labeled cell)
        CellwiseSourceParabolicPde<2> pde(cell_population, M_DUDT_COEFFICIENT,M_DIFFUSION_CONSTANT,M_UPTAKE_RATE);

        // Create a pde modifier and pass it to the simulation Add this first so in place for SimpleTargetArea one (calls cell pop update)
        MAKE_PTR_ARGS(ParabolicPDEModifier<2>, p_pde_modifier, (&pde));
        simulator.AddSimulationModifier(p_pde_modifier);

        // Add volume tracking modifier for CI cell cycle model.
        MAKE_PTR(VolumeTrackingModifier<2>, p_volume_tracking_modifier);
        simulator.AddSimulationModifier(p_volume_tracking_modifier);

        // A NagaiHondaForce has to be used together with an AbstractTargetAreaModifier
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);



        simulator.Solve();
    }

    void TestNodeBasedMonolayer() throw (Exception)
    {
        HoneycombMeshGenerator generator(M_NUM_CELLS_ACROSS, M_NUM_CELLS_ACROSS,0);
        MutableMesh<2,2>* p_generating_mesh = generator.GetMesh();
        NodesOnlyMesh<2>* p_mesh = new NodesOnlyMesh<2>;
        p_mesh->ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);

        std::vector<CellPtr> cells;
        GenerateCells(p_mesh->GetNumNodes(),cells,1.0);

        NodeBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("ParabolicMonolayers/Node");
        simulator.SetDt(1.0/120.0);
        simulator.SetSamplingTimestepMultiple(120);
        simulator.SetEndTime(M_TIME_FOR_SIMULATION);

        MAKE_PTR(RepulsionForce<2>, p_force);
        simulator.AddForce(p_force);

        // Set up PDE and pass to simulation via modifier (uniform secretion at each labeled cell)
        CellwiseSourceParabolicPde<2> pde(cell_population, M_DUDT_COEFFICIENT,M_DIFFUSION_CONSTANT,M_UPTAKE_RATE);

        // Create a pde modifier and pass it to the simulation
        MAKE_PTR_ARGS(ParabolicPDEModifier<2>, p_pde_modifier, (&pde));
        simulator.AddSimulationModifier(p_pde_modifier);

        // Add volume tracking modifier for CI cell cycle model.
        MAKE_PTR(VolumeTrackingModifier<2>, p_volume_tracking_modifier);
        simulator.AddSimulationModifier(p_volume_tracking_modifier);

        simulator.Solve();

        delete p_mesh; // to stop memory leaks
    }

    void TestMeshBasedMonolayer() throw (Exception)
    {
        HoneycombMeshGenerator generator(M_NUM_CELLS_ACROSS,M_NUM_CELLS_ACROSS,0);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        GenerateCells(p_mesh->GetNumNodes(),cells,1.0);

        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Set population to output all data to results files
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();

        cell_population.SetWriteVtkAsPoints(true);
        cell_population.AddPopulationWriter<VoronoiDataWriter>();


        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("ParabolicMonolayers/MeshPoint");
        simulator.SetDt(1.0/120.0);
        simulator.SetSamplingTimestepMultiple(120);
        simulator.SetEndTime(M_TIME_FOR_SIMULATION);

        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        p_force->SetCutOffLength(1.5);
        simulator.AddForce(p_force);

        // Set up PDE and pass to simulation via modifier (uniform secretion at each labeled cell)
        CellwiseSourceParabolicPde<2> pde(cell_population, M_DUDT_COEFFICIENT,M_DIFFUSION_CONSTANT,M_UPTAKE_RATE);

        // Create a pde modifier and pass it to the simulation
        MAKE_PTR_ARGS(ParabolicPDEModifier<2>, p_pde_modifier, (&pde));
        simulator.AddSimulationModifier(p_pde_modifier);

        // Add volume tracking modifier for CI cell cycle model.
        MAKE_PTR(VolumeTrackingModifier<2>, p_volume_tracking_modifier);
        simulator.AddSimulationModifier(p_volume_tracking_modifier);

        simulator.Solve();
    }

    void TestPottsBasedMonolayer() throw (Exception)
    {
        unsigned cell_width = 4;
        unsigned domain_width = 200;
        PottsMeshGenerator<2> generator(domain_width, M_NUM_CELLS_ACROSS, cell_width, domain_width, M_NUM_CELLS_ACROSS, cell_width);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;
        GenerateCells(p_mesh->GetNumElements(),cells,16);

        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetTemperature(0.1);
        cell_population.SetNumSweepsPerTimestep(100);

        // Set population to output all data to results files
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();

        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("ParabolicMonolayers/Potts");
        simulator.SetDt(1.0);
        simulator.SetSamplingTimestepMultiple(1);
        simulator.SetEndTime(M_TIME_FOR_SIMULATION);

        MAKE_PTR(VolumeConstraintPottsUpdateRule<2>, p_volume_constraint_update_rule);
        simulator.AddPottsUpdateRule(p_volume_constraint_update_rule);
        MAKE_PTR(SurfaceAreaConstraintPottsUpdateRule<2>, p_surface_area_update_rule);
        simulator.AddPottsUpdateRule(p_surface_area_update_rule);
        MAKE_PTR(AdhesionPottsUpdateRule<2>, p_adhesion_update_rule);
        simulator.AddPottsUpdateRule(p_adhesion_update_rule);

        // Set up PDE and pass to simulation via modifier (uniform secretion at each labeled cell)
        CellwiseSourceParabolicPde<2> pde(cell_population, M_DUDT_COEFFICIENT,M_DIFFUSION_CONSTANT,M_UPTAKE_RATE);

        // Create a pde modifier and pass it to the simulation
        MAKE_PTR_ARGS(ParabolicPDEModifier<2>, p_pde_modifier, (&pde));
        simulator.AddSimulationModifier(p_pde_modifier);

        // Add volume tracking modifier for CI cell cycle model.
        MAKE_PTR(VolumeTrackingModifier<2>, p_volume_tracking_modifier);
        simulator.AddSimulationModifier(p_volume_tracking_modifier);


        simulator.Solve();
    }

    void noTestCaBasedMonolayer() throw (Exception)
    {
        // Create a simple 2D PottsMesh
        unsigned domain_wide = 5*M_NUM_CELLS_ACROSS;

        PottsMeshGenerator<2> generator(domain_wide, 0, 0, domain_wide, 0, 0);
        PottsMesh<2>* p_mesh = generator.GetMesh();

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
        GenerateCells(location_indices.size(),cells,1);

        // Create cell population
        CaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);

        // Set population to output all data to results files
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellMutationStatesWriter>();

        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("ParabolicMonolayers/CA");
        simulator.SetDt(0.1);
        simulator.SetSamplingTimestepMultiple(10);
        simulator.SetEndTime(M_TIME_FOR_SIMULATION);

        // Adding update rule(s).
        MAKE_PTR(DiffusionCaUpdateRule<2u>, p_diffusion_update_rule);
        p_diffusion_update_rule->SetDiffusionParameter(0.1);

        simulator.AddCaUpdateRule(p_diffusion_update_rule);

        // Set up PDE and pass to simulation via modifier (uniform secretion at each labeled cell)
        CellwiseSourceParabolicPde<2> pde(cell_population, M_DUDT_COEFFICIENT,M_DIFFUSION_CONSTANT,M_UPTAKE_RATE);

        // Create a pde modifier and pass it to the simulation
        MAKE_PTR_ARGS(ParabolicPDEModifier<2>, p_pde_modifier, (&pde));
        simulator.AddSimulationModifier(p_pde_modifier);

       // Add volume tracking modifier for CI cell cycle model.
       MAKE_PTR(VolumeTrackingModifier<2>, p_volume_tracking_modifier);
       simulator.AddSimulationModifier(p_volume_tracking_modifier);

        // Run simulation
        simulator.Solve();
    }

};
