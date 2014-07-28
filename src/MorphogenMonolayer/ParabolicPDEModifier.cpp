/*

Copyright (c) 2005-2013, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "ParabolicPDEModifier.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "CaBasedCellPopulation.hpp"
#include "TetrahedralMesh.hpp"
#include "VtkMeshWriter.hpp"
#include "CellBasedParabolicPdeSolver.hpp"
#include "SimpleLinearParabolicSolver.hpp"
#include "AveragedSourcePde.hpp"
#include "Debug.hpp"

template<unsigned DIM>
ParabolicPDEModifier<DIM>::ParabolicPDEModifier(AbstractLinearParabolicPde<DIM,DIM>* pPde)
    : AbstractCellBasedSimulationModifier<DIM>(),
      mpPde(pPde),
      mSolution(NULL),
      mpFeMesh(NULL),
      mpBoundaryCondition(NULL),
      mDependentVariableName("morphogen"),
      mOutputDirectory("")
{
	assert(DIM==2);

	// Neumann BCS to be applied
	mpBoundaryCondition = new ConstBoundaryCondition<DIM>(0.0);
}

template<unsigned DIM>
ParabolicPDEModifier<DIM>::~ParabolicPDEModifier()
{
}


template<unsigned DIM>
void ParabolicPDEModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
	// Get the FeMesh from the cell population
	SetupFeMesh(rCellPopulation);

    // Set up boundary conditions
    std::auto_ptr<BoundaryConditionsContainer<DIM,DIM,1> > p_bcc = ConstructBoundaryConditionsContainer();

    // Construct the solution vector from cell data (takes care of cells dividing);
    UpdateSolutionVector(rCellPopulation);

    // Use CellBasedParabolicPdeSolver as cell wise PDE
    CellBasedParabolicPdeSolver<DIM> solver(mpFeMesh, mpPde, p_bcc.get());

    //TODO Investigate more than one pde timestep per spatial step
    SimulationTime* p_simulation_time = SimulationTime::Instance();
    double current_time = p_simulation_time->GetTime();
    double dt = p_simulation_time->GetTimeStep();
    solver.SetTimes(current_time,current_time + dt);
    solver.SetTimeStep(dt);

    // Use previous solution as ICS
    solver.SetInitialCondition(mSolution);

    mSolution = solver.Solve();

    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void ParabolicPDEModifier<DIM>::UpdateAtEndOfOutputTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    if (DIM>1)
    {
        std::ostringstream time_string;
        time_string << SimulationTime::Instance()->GetTimeStepsElapsed();
        std::string results_file = "pde_results_"+time_string.str();
        VtkMeshWriter<DIM,DIM>* p_vtk_mesh_writer = new VtkMeshWriter<DIM,DIM>(mOutputDirectory, results_file, false);

        ReplicatableVector solution_repl(mSolution);
        std::vector<double> pde_solution;
        for (unsigned i=0; i<mpFeMesh->GetNumNodes(); i++)
        {
           pde_solution.push_back(solution_repl[i]);
        }

        p_vtk_mesh_writer->AddPointData(mDependentVariableName,pde_solution);

        p_vtk_mesh_writer->WriteFilesUsingMesh(*mpFeMesh);
        delete p_vtk_mesh_writer;
    }
}

template<unsigned DIM>
void ParabolicPDEModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    // Cache the Output Directory
    mOutputDirectory = outputDirectory;

    //Setup a FeMesh to save the ics to
    SetupFeMesh(rCellPopulation);

    // Copy the cell data to mSolution (this is the initial conditions)
    UpdateSolutionVector(rCellPopulation);

    // Output the initial conditions on FeMesh
    UpdateAtEndOfOutputTimeStep(rCellPopulation);
}

template<unsigned DIM>
std::auto_ptr<BoundaryConditionsContainer<DIM,DIM,1> > ParabolicPDEModifier<DIM>::ConstructBoundaryConditionsContainer()
{
    std::auto_ptr<BoundaryConditionsContainer<DIM,DIM,1> > p_bcc(new BoundaryConditionsContainer<DIM,DIM,1>(false));

	for (typename TetrahedralMesh<DIM,DIM>::BoundaryElementIterator elem_iter = mpFeMesh->GetBoundaryElementIteratorBegin();
		 elem_iter != mpFeMesh->GetBoundaryElementIteratorEnd();
		 ++elem_iter)
	{
		p_bcc->AddNeumannBoundaryCondition(*elem_iter, mpBoundaryCondition);
	}

//	for (typename TetrahedralMesh<DIM,DIM>::BoundaryNodeIterator node_iter = mpFeMesh->GetBoundaryNodeIteratorBegin();
//	                 node_iter != mpFeMesh->GetBoundaryNodeIteratorEnd();
//	                 ++node_iter)
//	{
//		p_bcc->AddDirichletBoundaryCondition(*node_iter, mpBoundaryCondition);
//	}

	return p_bcc;
}


template<unsigned DIM>
void ParabolicPDEModifier<DIM>::SetupFeMesh(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
	// Make sure the cell population is in a nice state
	rCellPopulation.Update();

	// Get FE mesh from Cell Population
	if(dynamic_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation))
	{
		mpFeMesh = &(static_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation)->rGetMesh());
	}
	else if (dynamic_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation))
	{
		std::vector<Node<DIM> *> nodes;

		// Get the nodes of the NodesOnlyMesh
		for (typename AbstractMesh<DIM,DIM>::NodeIterator node_iter = rCellPopulation.rGetMesh().GetNodeIteratorBegin();
				 node_iter != rCellPopulation.rGetMesh().GetNodeIteratorEnd();
				 ++node_iter)
		{
				nodes.push_back(new Node<DIM>(node_iter->GetIndex(), node_iter->rGetLocation()));
		}

		mpFeMesh = new MutableMesh<DIM,DIM>(nodes);
		assert(mpFeMesh->GetNumNodes() == rCellPopulation.GetNumRealCells());

	}
	else if (dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation))
	{
		mpFeMesh = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation)->GetTetrahedralMeshUsingVertexMesh();
	}
	else if (dynamic_cast<PottsBasedCellPopulation<DIM>*>(&rCellPopulation))
	{
		std::vector<Node<DIM> *> nodes;

		// Create nodes at the centre of the cells
		for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
			 cell_iter != rCellPopulation.End();
			 ++cell_iter)
		{
			nodes.push_back(new Node<DIM>(rCellPopulation.GetLocationIndexUsingCell(*cell_iter), rCellPopulation.GetLocationOfCellCentre(*cell_iter)));
		}

		mpFeMesh = new MutableMesh<DIM,DIM>(nodes);
		assert(mpFeMesh->GetNumNodes() == rCellPopulation.GetNumRealCells());
	}
	else if (dynamic_cast<CaBasedCellPopulation<DIM>*>(&rCellPopulation))
	{
		std::vector<Node<DIM> *> nodes;

		// Create nodes at the centre of the cells
		unsigned cell_index = 0;
		for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
			 cell_iter != rCellPopulation.End();
			 ++cell_iter)
		{
			nodes.push_back(new Node<DIM>(cell_index, rCellPopulation.GetLocationOfCellCentre(*cell_iter)));
			cell_index++;
		}

		mpFeMesh = new MutableMesh<DIM,DIM>(nodes);
		assert(mpFeMesh->GetNumNodes() == rCellPopulation.GetNumRealCells());
	}
	else
	{
		NEVER_REACHED;
	}
}

template<unsigned DIM>
void ParabolicPDEModifier<DIM>::UpdateSolutionVector(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    //Clear and resize the solution vector
    PetscTools::Destroy(mSolution);
    mSolution = PetscTools::CreateAndSetVec(mpFeMesh->GetNumNodes(), 0.0);

    // Loop over nodes and get appropriate solution value from CellData


    //     Apply dirichlet Boundaries
    for (typename TetrahedralMesh<DIM,DIM>::NodeIterator node_iter = mpFeMesh->GetNodeIteratorBegin();
		  node_iter != mpFeMesh->GetNodeIteratorEnd();
		  ++node_iter)
    {
    	unsigned node_index = node_iter->GetIndex();
    	double solution_at_node;

    	if (dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation))
    	{
    		// Cells correspond to nodes in the Center of the vertex element
    		// nodes on vertices have averaged values from containing cells

    		unsigned num_vertex_nodes = rCellPopulation.GetNumNodes();
    		if (node_index >= num_vertex_nodes)
    		{
    			// Offset to relate elements in vertex mesh to nodes in tetrahedral mesh.

    			assert(node_index-num_vertex_nodes < num_vertex_nodes);
    			assert(node_index>=0);

    			CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(node_index-num_vertex_nodes);
    		    solution_at_node = p_cell->GetCellData()->GetItem(mDependentVariableName);
    		}
    		else
    		{
    			// Average over Data from containing elements (cells)
    			assert(node_index<num_vertex_nodes);
    			Node<DIM>* p_vertex_node = rCellPopulation.rGetMesh().GetNode(node_index);
    			std::set<unsigned> containing_elelments  = p_vertex_node->rGetContainingElementIndices();

    			solution_at_node = 0.0;

    			for (std::set<unsigned>::iterator index_iter = containing_elelments.begin();
		     		 index_iter != containing_elelments.end();
					 ++index_iter)
    			{

    				assert(*index_iter<num_vertex_nodes);
    				CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(*index_iter);
    				solution_at_node += p_cell->GetCellData()->GetItem(mDependentVariableName);
    			}
    			solution_at_node /= containing_elelments.size();
    		}
    	}
    	else
    	{
    		assert(dynamic_cast<AbstractCentreBasedCellPopulation<DIM>*>(&rCellPopulation)||
    			   dynamic_cast<AbstractOnLatticeCellPopulation<DIM>*>(&rCellPopulation));

			// Simple 1-1 correspondnece between cells and nodes in the FeMesh
			CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(node_index);
			solution_at_node = p_cell->GetCellData()->GetItem(mDependentVariableName);
		}

    	PetscVecTools::SetElement(mSolution,node_index,solution_at_node);
	}
}

template<unsigned DIM>
void ParabolicPDEModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Make sure the cell population is updated
    rCellPopulation.Update();

    // Store the PDE solution in an accessible form
    ReplicatableVector solution_repl(mSolution);

    // local cell index used by the CA simulation
    unsigned cell_index = 0;

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
                 cell_iter != rCellPopulation.End();
                 ++cell_iter)
	{
		unsigned tet_node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);

		if (dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation))
		{
			// Offset to relate elements in vertex mesh to nodes in tetrahedral mesh.
			tet_node_index += rCellPopulation.GetNumNodes();
		}

		if (dynamic_cast<CaBasedCellPopulation<DIM>*>(&rCellPopulation))
		{
			// here local cell index corresponds to tet node
			tet_node_index = cell_index;
			cell_index++;
		}

		double solution_at_node = solution_repl[tet_node_index];

		cell_iter->GetCellData()->SetItem(mDependentVariableName, solution_at_node);
	}
}

template<unsigned DIM>
void ParabolicPDEModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class ParabolicPDEModifier<1>;
template class ParabolicPDEModifier<2>;
template class ParabolicPDEModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ParabolicPDEModifier)

