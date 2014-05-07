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

#include "NeumannPDEModifier.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "MultipleCaBasedCellPopulation.hpp"
#include "TetrahedralMesh.hpp"
#include "VtkMeshWriter.hpp"
#include "CellBasedPdeSolver.hpp"
#include "SimpleLinearEllipticSolver.hpp"
#include "AveragedSourcePde.hpp"
#include "Debug.hpp"

template<unsigned DIM>
NeumannPDEModifier<DIM>::NeumannPDEModifier(AbstractLinearEllipticPde<DIM,DIM>* pPde)
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
NeumannPDEModifier<DIM>::~NeumannPDEModifier()
{
}


template<unsigned DIM>
void NeumannPDEModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
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
	else if (dynamic_cast<MultipleCaBasedCellPopulation<DIM>*>(&rCellPopulation))
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

    // Set up boundary conditions
    std::auto_ptr<BoundaryConditionsContainer<DIM,DIM,1> > p_bcc = ConstructBoundaryConditionsContainer();

    // Use CellBasedPdeSolver as cell wise PDE
    CellBasedPdeSolver<DIM> solver(mpFeMesh, mpPde, p_bcc.get());

    // TODO Use initial guess, when solving the system...
    mSolution = solver.Solve();

    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void NeumannPDEModifier<DIM>::UpdateAtEndOfOutputTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
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
void NeumannPDEModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    // Cache the Output Directory
    mOutputDirectory = outputDirectory;

    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    //UpdateCellData(rCellPopulation);

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
                    cell_iter != rCellPopulation.End();
                    ++cell_iter)
   	{
    	cell_iter->GetCellData()->SetItem(mDependentVariableName, 0.0);
   	}
}

template<unsigned DIM>
std::auto_ptr<BoundaryConditionsContainer<DIM,DIM,1> > NeumannPDEModifier<DIM>::ConstructBoundaryConditionsContainer()
{
    std::auto_ptr<BoundaryConditionsContainer<DIM,DIM,1> > p_bcc(new BoundaryConditionsContainer<DIM,DIM,1>(false));

//	for (typename TetrahedralMesh<DIM,DIM>::BoundaryElementIterator elem_iter = pMesh->GetBoundaryElementIteratorBegin();
//		 elem_iter != mpFeMesh->GetBoundaryElementIteratorEnd();
//		 ++elem_iter)
//	{
//		p_bcc->AddNeumannBoundaryCondition(*elem_iter, mpBoundaryCondition);
//	}
//
//	// Pin one node at zero so matrix is not singular
//	typename TetrahedralMesh<DIM,DIM>::BoundaryNodeIterator node_iter = mpFeMesh->GetBoundaryNodeIteratorBegin();
//	p_bcc->AddDirichletBoundaryCondition(*node_iter, mpBoundaryCondition);



	for (typename TetrahedralMesh<DIM,DIM>::BoundaryNodeIterator node_iter = mpFeMesh->GetBoundaryNodeIteratorBegin();
	                 node_iter != mpFeMesh->GetBoundaryNodeIteratorEnd();
	                 ++node_iter)
	{
		p_bcc->AddDirichletBoundaryCondition(*node_iter, mpBoundaryCondition);
	}

	return p_bcc;
}

template<unsigned DIM>
void NeumannPDEModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
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

		if (dynamic_cast<MultipleCaBasedCellPopulation<DIM>*>(&rCellPopulation))
		{
			// here local cell index corresponds to tet node
			tet_node_index = cell_index;
			cell_index++;
		}

		double solution_at_node = solution_repl[tet_node_index];

		cell_iter->GetCellData()->SetItem(mDependentVariableName, solution_at_node);
	}
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class NeumannPDEModifier<1>;
template class NeumannPDEModifier<2>;
template class NeumannPDEModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(NeumannPDEModifier)

