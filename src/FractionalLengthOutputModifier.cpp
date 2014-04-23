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

#include "FractionalLengthOutputModifier.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "Debug.hpp"

template<unsigned DIM>
FractionalLengthOutputModifier<DIM>::FractionalLengthOutputModifier()
    : AbstractCellBasedSimulationModifier<DIM>()
{
}

template<unsigned DIM>
FractionalLengthOutputModifier<DIM>::~FractionalLengthOutputModifier()
{
}

template<unsigned DIM>
void FractionalLengthOutputModifier<DIM>::UpdateAtEndOfOutputTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    CalculateFractionalLength(rCellPopulation);
}

template<unsigned DIM>
void FractionalLengthOutputModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    // Create output file
	OutputFileHandler output_file_handler( outputDirectory + "/", false);
    mpFractionalLengthsResultsFile = output_file_handler.OpenOutputFile("fractionallengths.dat");

    // Write Headers
    *mpFractionalLengthsResultsFile <<  "time \t fractional length \t total length \n";

    // Calculate before 1st timestep.
    CalculateFractionalLength(rCellPopulation);
}


template<unsigned DIM>
void FractionalLengthOutputModifier<DIM>::UpdateAtEndOfSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
	// Close output file.
    mpFractionalLengthsResultsFile->close();
}

template<unsigned DIM>
void FractionalLengthOutputModifier<DIM>::CalculateFractionalLength(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    double fractional_length = 0.0;
    double total_length = 0.0;
    double fractional_neighbours = 0.0;
    double total_neighbours = 0.0;

    // Make sure the cell population is updated
    rCellPopulation.Update();

    if (dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation))
    {
        MutableVertexMesh<DIM,DIM>* p_vertex_mesh = (&(dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation)->rGetMesh()));

        // Loop over vertex elements
        for (typename VertexMesh<DIM,DIM>::VertexElementIterator elem_iter = p_vertex_mesh->GetElementIteratorBegin();
             elem_iter != p_vertex_mesh->GetElementIteratorEnd();
             ++elem_iter)
        {
            unsigned element_index = elem_iter->GetIndex();

            // Get cell associated with this element
            CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(element_index);

            std::set<unsigned> neighbouring_elem_indices  = p_vertex_mesh->GetNeighbouringElementIndices(element_index);

            // Iterate over these neighbours
            for (std::set<unsigned>::iterator neighbour_iter = neighbouring_elem_indices.begin();
                    neighbour_iter != neighbouring_elem_indices.end();
                 ++neighbour_iter)
            {
                unsigned neighbour_index = *neighbour_iter;

                CellPtr p_neighbour_cell = rCellPopulation.GetCellUsingLocationIndex(neighbour_index);

                if ( (p_cell->template HasCellProperty<CellLabel>() && !(p_neighbour_cell->template HasCellProperty<CellLabel>())) ||
                     (!(p_cell->template HasCellProperty<CellLabel>()) && p_neighbour_cell->template HasCellProperty<CellLabel>()) )
                {
                    //Paranoia
                    assert((unsigned)(p_cell->template HasCellProperty<CellLabel>())+(unsigned)(p_neighbour_cell->template HasCellProperty<CellLabel>())==1);

                    fractional_length += p_vertex_mesh->GetEdgeLength(element_index,neighbour_index);
                    fractional_neighbours += 1.0;
                }

                total_length += p_vertex_mesh->GetEdgeLength(element_index,neighbour_index);
                total_neighbours += 1.0;
            }
        }
    }
    else if (dynamic_cast<PottsBasedCellPopulation<DIM>*>(&rCellPopulation))
	{
		PottsMesh<DIM>* p_potts_mesh = (&(dynamic_cast<PottsBasedCellPopulation<DIM>*>(&rCellPopulation)->rGetMesh()));

		// Loop over elements
		for (typename PottsMesh<DIM>::PottsElementIterator elem_iter = p_potts_mesh->GetElementIteratorBegin();
			 elem_iter != p_potts_mesh->GetElementIteratorEnd();
			 ++elem_iter)
		{
			unsigned element_index = elem_iter->GetIndex();

			// Get cell associated with this element
			CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(element_index);

			// loop over nodes in the element
			for (unsigned node_index = 0; node_index <  elem_iter->GetNumNodes(); node_index++)
			{
				unsigned node_global_index = elem_iter->GetNodeGlobalIndex(node_index);

				// Loop over neighbouring nodes Only want N,S,E,W neighbours as need to share an edge
				std::set<unsigned> neighbouring_node_indices  = p_potts_mesh->GetVonNeumannNeighbouringNodeIndices(node_global_index);

			   // Iterate over these neighbours
			   for (std::set<unsigned>::iterator neighbour_iter = neighbouring_node_indices.begin();
					   neighbour_iter != neighbouring_node_indices.end();
					++neighbour_iter)
			   {
				   std::set<unsigned> neighbouring_elem_indices = p_potts_mesh->GetNode(*neighbour_iter)->rGetContainingElementIndices();

				   assert(neighbouring_elem_indices.size()<2); // Either in element or in medium

				   if (neighbouring_elem_indices.size()==1)
				   {
					   unsigned neigbouring_elem_index =  *(neighbouring_elem_indices.begin());

					   if (neigbouring_elem_index != element_index)
					   {
						   // edge is between two different elements

						   CellPtr p_neighbour_cell = rCellPopulation.GetCellUsingLocationIndex(neigbouring_elem_index);

						   if ( (p_cell->template HasCellProperty<CellLabel>() && !(p_neighbour_cell->template HasCellProperty<CellLabel>())) ||
								(!(p_cell->template HasCellProperty<CellLabel>()) && p_neighbour_cell->template HasCellProperty<CellLabel>()) )
						   {
							   // Paranoia
							   assert((unsigned)(p_cell->template HasCellProperty<CellLabel>())+(unsigned)(p_neighbour_cell->template HasCellProperty<CellLabel>())==1);

							   fractional_length += 1.0;
						   }
						   total_length += 1.0;
					   }
				   }
				   else
				   {
					   // Original node is on boundary of mesh so edge only counted once
					   total_length += 2.0;
				   }
			   }

			}
		}
	}
    else if (dynamic_cast<MeshBasedCellPopulationWithGhostNodes<DIM>*>(&rCellPopulation))
    {
    	MeshBasedCellPopulationWithGhostNodes<DIM>* p_population = static_cast<MeshBasedCellPopulationWithGhostNodes<DIM>*>(&rCellPopulation);

    	p_population->CreateVoronoiTessellation();

        // Loop over Cells
		for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
			 cell_iter != rCellPopulation.End();
			 ++cell_iter)
		{
			// Get the location index corresponding to this cell
			unsigned index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);

			// Get the set of neighbouring location indices
			std::set<unsigned> neighbour_indices = rCellPopulation.GetNeighbouringNodeIndices(index);

            // Iterate over these neighbours
            for (std::set<unsigned>::iterator neighbour_iter = neighbour_indices.begin();
                    neighbour_iter != neighbour_indices.end();
                 ++neighbour_iter)
            {
                unsigned neighbour_index = *neighbour_iter;

                //PRINT_2_VARIABLES(index,neighbour_index);
                double edge_length = p_population->GetVoronoiEdgeLength(index,neighbour_index);

				if ( p_population->IsGhostNode(neighbour_index) )
				{
					// On the edge so account for edge twice
					total_length +=2.0*edge_length;
				}
				else
				{
					total_length +=edge_length;
					total_neighbours +=1;

					CellPtr p_neighbour_cell = rCellPopulation.GetCellUsingLocationIndex(neighbour_index);

					if ( (cell_iter->template HasCellProperty<CellLabel>() && !(p_neighbour_cell->template HasCellProperty<CellLabel>())) ||
						 (!(cell_iter->template HasCellProperty<CellLabel>()) && p_neighbour_cell->template HasCellProperty<CellLabel>()) )
					{
						fractional_length += edge_length;
						fractional_neighbours += 1.0;
					}
				}
            }
		}
    }
    else if (dynamic_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation))
	{
		// Loop over Cells
		for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
			 cell_iter != rCellPopulation.End();
			 ++cell_iter)
		{
			// Get the location index corresponding to this cell
			unsigned index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);

			// Get the set of neighbouring location indices
			std::set<unsigned> neighbour_indices;
			assert(dynamic_cast<AbstractCentreBasedCellPopulation<DIM>*>(&rCellPopulation));
			{
				neighbour_indices = rCellPopulation.GetNeighbouringNodeIndices(index);
			}

			if (!neighbour_indices.empty())
			{
				for (std::set<unsigned>::iterator iter = neighbour_indices.begin();
					 iter != neighbour_indices.end();
					 ++iter)
				{
					CellPtr p_neighbour_cell = rCellPopulation.GetCellUsingLocationIndex(*iter);

					//Approximate contact area between cells
					double cell_seperation = norm_2(rCellPopulation.GetLocationOfCellCentre(*cell_iter) - rCellPopulation.GetLocationOfCellCentre(p_neighbour_cell));
					double effective_cell_radius = 0.6;
					double contact_area = 0.0;

					if 	(cell_seperation < 2.0*	effective_cell_radius)
					{
						contact_area = 2.0*sqrt(effective_cell_radius*effective_cell_radius - cell_seperation*cell_seperation/4.0 );
					}
					// else not in contact so zero

					total_length +=contact_area;
					total_neighbours +=1;

					if ( (cell_iter->template HasCellProperty<CellLabel>() && !(p_neighbour_cell->template HasCellProperty<CellLabel>())) ||
						 (!(cell_iter->template HasCellProperty<CellLabel>()) && p_neighbour_cell->template HasCellProperty<CellLabel>()) )
					{
						//Paranoia
						assert((unsigned)(cell_iter->template HasCellProperty<CellLabel>())+(unsigned)(p_neighbour_cell->template HasCellProperty<CellLabel>())==1);
						fractional_length += contact_area;
						fractional_neighbours += 1.0;
					}
				}

			}
		}
	}

    else
    {
        EXCEPTION("FractionalLengthOffLatticeSimulation only works for Vertex Potts and Node based simulations at present");
    }

    // each edge is counted twice so divide by 2
    fractional_length /= 2.0;
    total_length /= 2.0;

    // each neighbour is counted twice so divide by 2
	fractional_neighbours /= 2.0;
	total_neighbours /= 2.0;




    *mpFractionalLengthsResultsFile <<  SimulationTime::Instance()->GetTime() << "\t" << fractional_length << "\t" << total_length << "\t" << fractional_neighbours << "\t" << total_neighbours <<"\n";
}



/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class FractionalLengthOutputModifier<1>;
template class FractionalLengthOutputModifier<2>;
template class FractionalLengthOutputModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(FractionalLengthOutputModifier)

