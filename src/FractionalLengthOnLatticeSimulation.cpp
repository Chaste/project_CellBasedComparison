/*

Copyright (c) 2005-2012, University of Oxford.
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

#include "FractionalLengthOnLatticeSimulation.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "Debug.hpp"

template<unsigned DIM>
FractionalLengthOnLatticeSimulation<DIM>::FractionalLengthOnLatticeSimulation(AbstractCellPopulation<DIM>& rCellPopulation,
                                                                  bool deleteCellPopulationInDestructor,
                                                                  bool initialiseCells)
    : OnLatticeSimulation<DIM>(rCellPopulation, deleteCellPopulationInDestructor, initialiseCells)
{
}

template<unsigned DIM>
FractionalLengthOnLatticeSimulation<DIM>::~FractionalLengthOnLatticeSimulation()
{
}

template<unsigned DIM>
void FractionalLengthOnLatticeSimulation<DIM>::SetupSolve()
{
    // First call method on base class
    OnLatticeSimulation<DIM>::SetupSolve();

    // Make Output File
    OutputFileHandler output_file_handler(this->mSimulationOutputDirectory + "/", false);
    mFractionalLengthsResultsFile = output_file_handler.OpenOutputFile("fractionallengths.dat");

    CalculateFractionalLength();
}

template<unsigned DIM>
void FractionalLengthOnLatticeSimulation<DIM>::UpdateAtEndOfSolve()
{
    // First call method on base class
    OnLatticeSimulation<DIM>::UpdateAtEndOfSolve();

    mFractionalLengthsResultsFile->close();
}

template<unsigned DIM>
void FractionalLengthOnLatticeSimulation<DIM>::UpdateAtEndOfTimeStep()
{
    // First call method on base class
    OnLatticeSimulation<DIM>::UpdateAtEndOfTimeStep();

    CalculateFractionalLength();
}

template<unsigned DIM>
void FractionalLengthOnLatticeSimulation<DIM>::CalculateFractionalLength()
{
    double fractional_length = 0.0;
    double total_length = 0.0;

    // Make sure the cell population is updated
    this->mrCellPopulation.Update();

    if (dynamic_cast<PottsBasedCellPopulation<DIM>*>(&this->mrCellPopulation))
    {
        PottsMesh<DIM>* p_potts_mesh = (&(dynamic_cast<PottsBasedCellPopulation<DIM>*>(&(this->mrCellPopulation))->rGetMesh()));

        // Loop over elements
        for (typename PottsMesh<DIM>::PottsElementIterator elem_iter = p_potts_mesh->GetElementIteratorBegin();
             elem_iter != p_potts_mesh->GetElementIteratorEnd();
             ++elem_iter)
        {
            unsigned element_index = elem_iter->GetIndex();

            // Get cell associated with this element
            CellPtr p_cell = this->mrCellPopulation.GetCellUsingLocationIndex(element_index);

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

                           CellPtr p_neighbour_cell = this->mrCellPopulation.GetCellUsingLocationIndex(neigbouring_elem_index);

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
    else
    {
        EXCEPTION("FractionalLengthOnLatticeSimulation only works for Vertex simulations at present");
    }

    // each edge is counted twice so divide by 2
    fractional_length /= 2.0;
    total_length /= 2.0;

    *mFractionalLengthsResultsFile <<  SimulationTime::Instance()->GetTime() << "\t" << fractional_length << "\t" << total_length <<"\n";
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class FractionalLengthOnLatticeSimulation<1>;
template class FractionalLengthOnLatticeSimulation<2>;
template class FractionalLengthOnLatticeSimulation<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(FractionalLengthOnLatticeSimulation)