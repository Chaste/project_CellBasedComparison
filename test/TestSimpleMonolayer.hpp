#include <cxxtest/TestSuite.h>
#include "CellBasedSimulationArchiver.hpp"
#include "SmartPointers.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"

#include "TransitCellProliferativeType.hpp"

#include "StochasticDurationCellCycleModel.hpp"
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
#include "ShovingCaBasedDivisionRule.hpp"
#include "CryptShovingCaBasedDivisionRule.hpp"

#include "CellAgesWriter.hpp"
#include "CellIdWriter.hpp"

#include "PetscSetupAndFinalize.hpp"

/*
 * Just simulates a free growing monolayer.
 */
class TestSimpleMonolayers : public AbstractCellBasedWithTimingsTestSuite
{
public:
    void NoTestVertexBasedMonolayer() throw (Exception)
    {
        // Create Mesh
        HoneycombVertexMeshGenerator generator(2, 2);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();

        MAKE_PTR(TransitCellProliferativeType, p_transit_type);

        // Create Cells
        std::vector<CellPtr> cells;
        CellsGenerator<StochasticDurationCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_transit_type);

        // Create Population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // Create Simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("Monolayers/Vertex");
        simulator.SetSamplingTimestepMultiple(200);
        simulator.SetEndTime(20.0);

        // Create Forces and pass to simulation
        MAKE_PTR(NagaiHondaForce<2>, p_force);
        simulator.AddForce(p_force);

        // A NagaiHondaForce has to be used together with an AbstractTargetAreaModifier
        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        // Create CellKillers and pass to simulation

        // Create CellPopulationBoundaryCondiditons and pass to simulation

        simulator.Solve();
    }

    void NoTestNodeBasedMonolayer() throw (Exception)
    {
        HoneycombMeshGenerator generator(2, 2);
        MutableMesh<2,2>* p_generating_mesh = generator.GetMesh();
        NodesOnlyMesh<2>* p_mesh = new NodesOnlyMesh<2>;
        p_mesh->ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);


        MAKE_PTR(TransitCellProliferativeType, p_transit_type);

        std::vector<CellPtr> cells;
        CellsGenerator<StochasticDurationCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), p_transit_type);

        NodeBasedCellPopulation<2> cell_population(*p_mesh, cells);

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("Monolayers/Node");
        simulator.SetSamplingTimestepMultiple(12);
        simulator.SetEndTime(20.0);

        MAKE_PTR(RepulsionForce<2>, p_force);
        simulator.AddForce(p_force);

        simulator.Solve();

        delete p_mesh; // to stop memory leaks
    }

    void NoTestMeshBasedMonolayer() throw (Exception)
    {
        HoneycombMeshGenerator generator(2, 2);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();

        MAKE_PTR(TransitCellProliferativeType, p_transit_type);

        std::vector<CellPtr> cells;
        CellsGenerator<StochasticDurationCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), p_transit_type);

        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("Monolayers/Mesh");
        simulator.SetSamplingTimestepMultiple(12);
        simulator.SetEndTime(20.0);

        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        simulator.AddForce(p_force);

        simulator.Solve();
    }

    void NoTestPottsBasedMonolayer() throw (Exception)
    {
        PottsMeshGenerator<2> generator(20, 2, 4, 20, 2, 4);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        MAKE_PTR(TransitCellProliferativeType, p_transit_type);

        std::vector<CellPtr> cells;
        CellsGenerator<StochasticDurationCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_transit_type);

        PottsBasedCellPopulation<2> cell_population(*p_mesh, cells);
        cell_population.SetTemperature(1.0);

        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("Monolayers/Potts");
        simulator.SetEndTime(20.0);

        MAKE_PTR(VolumeConstraintPottsUpdateRule<2>, p_volume_constraint_update_rule);
        simulator.AddPottsUpdateRule(p_volume_constraint_update_rule);
        MAKE_PTR(SurfaceAreaConstraintPottsUpdateRule<2>, p_surface_area_update_rule);
        simulator.AddPottsUpdateRule(p_surface_area_update_rule);
        MAKE_PTR(AdhesionPottsUpdateRule<2>, p_adhesion_update_rule);
        simulator.AddPottsUpdateRule(p_adhesion_update_rule);

        simulator.Solve();
    }

    void TestCaBasedMonolayer() throw (Exception)
    {
        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(100, 0, 0, 100, 0, 0);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Specify where cells lie
        std::vector<unsigned> location_indices;
        location_indices.push_back(5050u);

        MAKE_PTR(TransitCellProliferativeType, p_transit_type);

        std::vector<CellPtr> cells;
        CellsGenerator<StochasticDurationCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, location_indices.size(), p_transit_type);

        // Create cell population
        CaBasedCellPopulation<2> cell_population(*p_mesh, cells, location_indices);

        boost::shared_ptr<AbstractCaBasedDivisionRule<2> > p_division_rule(new ShovingCaBasedDivisionRule<2>());
        cell_population.SetCaBasedDivisionRule(p_division_rule);



        cell_population.AddCellWriter<CellAgesWriter>();
        cell_population.AddCellWriter<CellIdWriter>();

        OnLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("Monolayers/CA");
        simulator.SetDt(1.0);
        simulator.SetEndTime(150);

        // Adding update rule(s).
//        MAKE_PTR(DiffusionCaUpdateRule<2u>, p_diffusion_update_rule);
//        p_diffusion_update_rule->SetDiffusionParameter(0.05);
//        simulator.AddCaUpdateRule(p_diffusion_update_rule);

        // Run simulation
        simulator.Solve();
    }

};
