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

#include "AdjacencyMatrixOutputModifier.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "Debug.hpp"

template<unsigned DIM>
AdjacencyMatrixOutputModifier<DIM>::AdjacencyMatrixOutputModifier()
    : AbstractCellBasedSimulationModifier<DIM>()
{
}

template<unsigned DIM>
AdjacencyMatrixOutputModifier<DIM>::~AdjacencyMatrixOutputModifier()
{
}

template<unsigned DIM>
void AdjacencyMatrixOutputModifier<DIM>::UpdateAtEndOfOutputTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    CalculateAdjacencyMatix(rCellPopulation);
}

template<unsigned DIM>
void AdjacencyMatrixOutputModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    // Create output file
    OutputFileHandler output_file_handler( outputDirectory + "/", false);
    mpAdjacencyMatrixResultsFile = output_file_handler.OpenOutputFile("adjacencymatrices.dat");

    // Write Headers
    //*mpAdjacencyMatrixResultsFile <<  "time \t number of cells \t adjacency matrix \n";

    // Calculate before 1st timestep.
    CalculateAdjacencyMatix(rCellPopulation);
}


template<unsigned DIM>
void AdjacencyMatrixOutputModifier<DIM>::UpdateAtEndOfSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Close output file.
    mpAdjacencyMatrixResultsFile->close();
}

template<unsigned DIM>
void AdjacencyMatrixOutputModifier<DIM>::CalculateAdjacencyMatix(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Make sure the cell population is updated
    rCellPopulation.Update();

    unsigned num_cells = rCellPopulation.GetNumRealCells();

    // Store a vector to map between cells numbered 1-n and locations indices
    std::map<unsigned,unsigned> local_cell_id_location_index_map;

    unsigned local_cell_id = 0;
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End();
             ++cell_iter)
        {
            local_cell_id_location_index_map[rCellPopulation.GetLocationIndexUsingCell(*cell_iter)] = local_cell_id;
            local_cell_id++;
        }
    assert(local_cell_id = num_cells+1);

    // Loop over all cells and calculate the adjacency matrix
    unsigned adjacency_matrix[num_cells*num_cells];
    for (unsigned i=0; i<num_cells*num_cells; i++)
    {
        adjacency_matrix[i] = 0;
    }

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        // Get the location index corresponding to this cell
        unsigned index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);

        // Get the set of neighbouring location indices
        std::set<unsigned> neighbour_indices;

        if (dynamic_cast<MeshBasedCellPopulationWithGhostNodes<DIM>*>(&rCellPopulation))
        {
            MeshBasedCellPopulationWithGhostNodes<DIM> * p_cell_population = static_cast<MeshBasedCellPopulationWithGhostNodes<DIM>*>(&rCellPopulation);

            neighbour_indices = rCellPopulation.GetNeighbouringNodeIndices(index);

            // Remove ghost nodes from the neighbour indices
            for(std::set<unsigned>::iterator iter = neighbour_indices.begin();
                iter != neighbour_indices.end();
                )
            {
               if (p_cell_population->IsGhostNode(*iter))
               {
                   neighbour_indices.erase(iter++);
               }
               else
               {
                  ++iter;
               }
            }
        }
        else if (dynamic_cast<AbstractCentreBasedCellPopulation<DIM>*>(&rCellPopulation))
        {
            neighbour_indices = rCellPopulation.GetNeighbouringNodeIndices(index);
        }
        else if (dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation))
        {
            neighbour_indices = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation)->rGetMesh().GetNeighbouringElementIndices(index);
        }
        else if(dynamic_cast<PottsBasedCellPopulation<DIM>*>(&rCellPopulation))
        {
            neighbour_indices = static_cast<PottsBasedCellPopulation<DIM>*>(&rCellPopulation)->rGetMesh().GetNeighbouringElementIndices(index);
        }
        else
        {
            ///\todo #2265 - this functionality is not yet implemented in ca-based simulations
            NEVER_REACHED;
        }

        if (!neighbour_indices.empty())
        {
            unsigned local_cell_index = local_cell_id_location_index_map[index];

            for (std::set<unsigned>::iterator neighbour_iter = neighbour_indices.begin();
                 neighbour_iter != neighbour_indices.end();
                 ++neighbour_iter)
            {
                // 1 - non non; 2 -label label; 3 - label non
                unsigned type_of_link = 1;

                CellPtr p_neighbour_cell = rCellPopulation.GetCellUsingLocationIndex(*neighbour_iter);

                if ( (cell_iter->template HasCellProperty<CellLabel>() && !(p_neighbour_cell->template HasCellProperty<CellLabel>())) ||
                     (!(cell_iter->template HasCellProperty<CellLabel>()) && p_neighbour_cell->template HasCellProperty<CellLabel>()) )
                {
                    type_of_link = 3;
                }
                else if (cell_iter->template HasCellProperty<CellLabel>() && p_neighbour_cell->template HasCellProperty<CellLabel>())
                {
                    type_of_link = 2;
                }
                else
                {
                    assert (!(cell_iter->template HasCellProperty<CellLabel>()) && !(p_neighbour_cell->template HasCellProperty<CellLabel>()));
                }

                unsigned local_neighbour_index = local_cell_id_location_index_map[*neighbour_iter];
                adjacency_matrix[local_cell_index + num_cells*local_neighbour_index] = type_of_link;
                adjacency_matrix[num_cells*local_cell_index + local_neighbour_index] = type_of_link;
            }
        }
        else
        {
            // If simulation trips at this assertion it is because at least one of the cells has no neighbours (as defined by mesh/population/interaction distance)
            //NEVER_REACHED;
        }
    }

    *mpAdjacencyMatrixResultsFile <<  SimulationTime::Instance()->GetTime() << "\t" << num_cells << "\t";
    for (unsigned i=0; i<num_cells*num_cells; i++)
    {
        *mpAdjacencyMatrixResultsFile << adjacency_matrix[i] << "\t";
    }
    *mpAdjacencyMatrixResultsFile << "\n";
}

// Explicit instantiation
template class AdjacencyMatrixOutputModifier<1>;
template class AdjacencyMatrixOutputModifier<2>;
template class AdjacencyMatrixOutputModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(AdjacencyMatrixOutputModifier)

