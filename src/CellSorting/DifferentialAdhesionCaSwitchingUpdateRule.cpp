/*

Copyright (c) 2005-2015, University of Oxford.
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

#include "DifferentialAdhesionCaSwitchingUpdateRule.hpp"

template<unsigned DIM>
DifferentialAdhesionCaSwitchingUpdateRule<DIM>::DifferentialAdhesionCaSwitchingUpdateRule()
    : AbstractCaSwitchingUpdateRule<DIM>(),
      mCellCellAdhesionEnergyParameter(0.2), // Educated guess
      mLabelledCellLabelledCellAdhesionEnergyParameter(0.2), // Educated guess
      mLabelledCellCellAdhesionEnergyParameter(0.1), // Educated guess
      mCellBoundaryAdhesionEnergyParameter(0.2), // Educated guess
      mLabelledCellBoundaryAdhesionEnergyParameter(0.2) // Educated guess
{
}

template<unsigned DIM>
DifferentialAdhesionCaSwitchingUpdateRule<DIM>::~DifferentialAdhesionCaSwitchingUpdateRule()
{
}

template<unsigned DIM>
double DifferentialAdhesionCaSwitchingUpdateRule<DIM>::EvaluateSwitchingProbability(unsigned currentNodeIndex,
                                                                      unsigned neighbourNodeIndex,
                                                                      CaBasedCellPopulation<DIM>& rCellPopulation,
                                                                      double dt,
                                                                      double deltaX)
{
    CellPtr p_cell_1 = rCellPopulation.GetCellUsingLocationIndex(currentNodeIndex);
    CellPtr p_cell_2 = rCellPopulation.GetCellUsingLocationIndex(neighbourNodeIndex);
    bool is_cell_1_labelled = p_cell_1->HasCellProperty<CellLabel>();
    bool is_cell_2_labelled = p_cell_2->HasCellProperty<CellLabel>();

    // Energy before and after switch
    double H_0=0.0;
    double H_1=0.0;

    std::set<unsigned> node_1_neighbouring_node_indices = rCellPopulation.rGetMesh().GetMooreNeighbouringNodeIndices(currentNodeIndex);

    for (std::set<unsigned>::iterator iter = node_1_neighbouring_node_indices.begin();
                 iter != node_1_neighbouring_node_indices.end();
                 ++iter)
    {
        if(rCellPopulation.IsCellAttachedToLocationIndex(*iter))
        {
            // Cell attached to neighbour cell

            bool is_neighbour_labelled = rCellPopulation.GetCellUsingLocationIndex(*iter)->template HasCellProperty<CellLabel>();

            if (is_neighbour_labelled)
            {
                if (is_cell_1_labelled)
                {
                    H_0 += mLabelledCellLabelledCellAdhesionEnergyParameter;
                }
                else
                {
                    H_0 += mLabelledCellCellAdhesionEnergyParameter;
                }
                if (is_cell_2_labelled)
                {
                    H_1 += mLabelledCellLabelledCellAdhesionEnergyParameter;
                }
                else
                {
                    H_1 += mLabelledCellCellAdhesionEnergyParameter;
                }

            }
            else
            {
                if (is_cell_1_labelled)
                {
                    H_0 += mLabelledCellCellAdhesionEnergyParameter;
                }
                else
                {
                    H_0 += mCellCellAdhesionEnergyParameter;
                }
                if (is_cell_2_labelled)
                {
                    H_1 += mLabelledCellCellAdhesionEnergyParameter;
                }
                else
                {
                    H_1 += mCellCellAdhesionEnergyParameter;
                }
            }
        }
        else // No cell on neighbour
        {
            if (is_cell_1_labelled)
            {
                H_0 += mLabelledCellBoundaryAdhesionEnergyParameter;
            }
            else
            {
                H_0 += mCellBoundaryAdhesionEnergyParameter;
            }
            if (is_cell_2_labelled)
            {
                H_1 += mLabelledCellBoundaryAdhesionEnergyParameter;
            }
            else
            {
                H_1 += mCellBoundaryAdhesionEnergyParameter;
            }
        }
    }


    std::set<unsigned> node_2_neighbouring_node_indices = rCellPopulation.rGetMesh().GetMooreNeighbouringNodeIndices(neighbourNodeIndex);

    for (std::set<unsigned>::iterator iter = node_2_neighbouring_node_indices.begin();
                 iter != node_2_neighbouring_node_indices.end();
                 ++iter)
    {
        if(rCellPopulation.IsCellAttachedToLocationIndex(*iter))
        {
            // Cell attached to neighbour cell

            bool is_neighbour_labelled = rCellPopulation.GetCellUsingLocationIndex(*iter)->template HasCellProperty<CellLabel>();

            if (is_neighbour_labelled)
            {
                if (is_cell_1_labelled)
                {
                    H_1 += mLabelledCellLabelledCellAdhesionEnergyParameter;
                }
                else
                {
                    H_1 += mLabelledCellCellAdhesionEnergyParameter;
                }
                if (is_cell_2_labelled)
                {
                    H_0 += mLabelledCellLabelledCellAdhesionEnergyParameter;
                }
                else
                {
                    H_0 += mLabelledCellCellAdhesionEnergyParameter;
                }

            }
            else
            {
                if (is_cell_1_labelled)
                {
                    H_1 += mLabelledCellCellAdhesionEnergyParameter;
                }
                else
                {
                    H_1 += mCellCellAdhesionEnergyParameter;
                }
                if (is_cell_2_labelled)
                {
                    H_0 += mLabelledCellCellAdhesionEnergyParameter;
                }
                else
                {
                    H_0 += mCellCellAdhesionEnergyParameter;
                }
            }
        }
        else // No cell on neighbour
        {
            if (is_cell_1_labelled)
            {
                H_1 += mLabelledCellBoundaryAdhesionEnergyParameter;
            }
            else
            {
                H_1 += mCellBoundaryAdhesionEnergyParameter;
            }
            if (is_cell_2_labelled)
            {
                H_0 += mLabelledCellBoundaryAdhesionEnergyParameter;
            }
            else
            {
                H_0 += mCellBoundaryAdhesionEnergyParameter;
            }
        }
    }

    double probability_of_switch = 0.0;

    if (H_1-H_0>0)
    {
        probability_of_switch =  (H_1-H_0)*dt;
    }

     return probability_of_switch;


}

template<unsigned DIM>
double DifferentialAdhesionCaSwitchingUpdateRule<DIM>::GetCellCellAdhesionEnergyParameter()
{
    return mCellCellAdhesionEnergyParameter;
}

template<unsigned DIM>
void DifferentialAdhesionCaSwitchingUpdateRule<DIM>::SetCellCellAdhesionEnergyParameter(double cellCellAdhesionEnergyParameter)
{
    mCellCellAdhesionEnergyParameter = cellCellAdhesionEnergyParameter;
}

template<unsigned DIM>
double DifferentialAdhesionCaSwitchingUpdateRule<DIM>::GetLabelledCellLabelledCellAdhesionEnergyParameter()
{
    return mLabelledCellLabelledCellAdhesionEnergyParameter;
}

template<unsigned DIM>
void DifferentialAdhesionCaSwitchingUpdateRule<DIM>::SetLabelledCellLabelledCellAdhesionEnergyParameter(double labelledCellLabelledCellAdhesionEnergyParameter)
{
    mLabelledCellLabelledCellAdhesionEnergyParameter = labelledCellLabelledCellAdhesionEnergyParameter;
}

template<unsigned DIM>
double DifferentialAdhesionCaSwitchingUpdateRule<DIM>::GetLabelledCellCellAdhesionEnergyParameter()
{
    return mLabelledCellCellAdhesionEnergyParameter;
}

template<unsigned DIM>
void DifferentialAdhesionCaSwitchingUpdateRule<DIM>::SetLabelledCellCellAdhesionEnergyParameter(double labelledCellCellAdhesionEnergyParameter)
{
    mLabelledCellCellAdhesionEnergyParameter = labelledCellCellAdhesionEnergyParameter;
}

template<unsigned DIM>
double DifferentialAdhesionCaSwitchingUpdateRule<DIM>::GetCellBoundaryAdhesionEnergyParameter()
{
    return mCellBoundaryAdhesionEnergyParameter;
}

template<unsigned DIM>
void DifferentialAdhesionCaSwitchingUpdateRule<DIM>::SetCellBoundaryAdhesionEnergyParameter(double cellBoundaryAdhesionEnergyParameter)
{
    mCellBoundaryAdhesionEnergyParameter = cellBoundaryAdhesionEnergyParameter;
}

template<unsigned DIM>
double DifferentialAdhesionCaSwitchingUpdateRule<DIM>::GetLabelledCellBoundaryAdhesionEnergyParameter()
{
    return mLabelledCellBoundaryAdhesionEnergyParameter;
}

template<unsigned DIM>
void DifferentialAdhesionCaSwitchingUpdateRule<DIM>::SetLabelledCellBoundaryAdhesionEnergyParameter(double labelledCellBoundaryAdhesionEnergyParameter)
{
    mLabelledCellBoundaryAdhesionEnergyParameter = labelledCellBoundaryAdhesionEnergyParameter;
}

template<unsigned DIM>
void DifferentialAdhesionCaSwitchingUpdateRule<DIM>::OutputSwitchingUpdateRuleParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<CellCellAdhesionEnergyParameter>" << mCellCellAdhesionEnergyParameter << "</CellCellAdhesionEnergyParameter>\n";
    *rParamsFile << "\t\t\t<LabelledCellLabelledCellAdhesionEnergyParameter>" << mLabelledCellLabelledCellAdhesionEnergyParameter << "</LabelledCellLabelledCellAdhesionEnergyParameter>\n";
    *rParamsFile << "\t\t\t<LabelledCellCellAdhesionEnergyParameter>" << mLabelledCellCellAdhesionEnergyParameter << "</LabelledCellCellAdhesionEnergyParameter>\n";
    *rParamsFile << "\t\t\t<CellBoundaryAdhesionEnergyParameter>" << mCellBoundaryAdhesionEnergyParameter << "</CellBoundaryAdhesionEnergyParameter>\n";
    *rParamsFile << "\t\t\t<LabelledCellBoundaryAdhesionEnergyParameter>" << mLabelledCellBoundaryAdhesionEnergyParameter << "</LabelledCellBoundaryAdhesionEnergyParameter>\n";

    // Call method on direct parent class
    AbstractCaSwitchingUpdateRule<DIM>::OutputSwitchingUpdateRuleParameters(rParamsFile);
}

// Explicit instantiation
template class DifferentialAdhesionCaSwitchingUpdateRule<1u>;
template class DifferentialAdhesionCaSwitchingUpdateRule<2u>;
template class DifferentialAdhesionCaSwitchingUpdateRule<3u>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(DifferentialAdhesionCaSwitchingUpdateRule)
