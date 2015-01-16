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

#include "StemCellRetainerForce.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "StemCellProliferativeType.hpp"
#include "PanethCellProliferativeType.hpp"

template<unsigned DIM>
StemCellRetainerForce<DIM>::StemCellRetainerForce()
    : AbstractForce<DIM>(),
      mStemCellForceMagnitudeParameter(10.0)
{
}

template<unsigned DIM>
StemCellRetainerForce<DIM>::~StemCellRetainerForce()
{
}

template<unsigned DIM>
void StemCellRetainerForce<DIM>::SetStemCellForceMagnitudeParameter(double stemCellForceMagnitudeParameter)
{
    mStemCellForceMagnitudeParameter = stemCellForceMagnitudeParameter;
}

template<unsigned DIM>
double StemCellRetainerForce<DIM>::GetStemCellForceMagnitudeParameter()
{
    return mStemCellForceMagnitudeParameter;
}

template<unsigned DIM>
void StemCellRetainerForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{

    if (dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation))
    {
        for (unsigned node_index=0; node_index<rCellPopulation.GetNumNodes(); node_index++)
        {
            Node<DIM>* p_node = rCellPopulation.GetNode(node_index);

            bool node_in_stem_cell  = false;

            std::set<unsigned> containing_element_indices = p_node->rGetContainingElementIndices();

            for (std::set<unsigned>::iterator element_iter = containing_element_indices.begin();
                    element_iter != containing_element_indices.end();
                    ++element_iter)
            {
                 boost::shared_ptr<AbstractCellProliferativeType> p_cell_type = rCellPopulation.GetCellUsingLocationIndex(*element_iter)->GetCellProliferativeType();


                if (p_cell_type->IsType<StemCellProliferativeType>())
                {
                    node_in_stem_cell = true;
                    break;
                }
            }

            if (node_in_stem_cell)
            {
                double height = p_node->rGetLocation()[DIM-1];

                c_vector<double,DIM> force_direction = zero_vector<double>(DIM);
                force_direction[DIM-1]=-1.0;

                c_vector<double, DIM> force_contribution =mStemCellForceMagnitudeParameter*force_direction*(1.0+height*height);
                p_node->AddAppliedForceContribution(force_contribution);
            }
        }
    }
    else // A centre based simulation
    {
        assert(dynamic_cast<AbstractCentreBasedCellPopulation<DIM>*>(&rCellPopulation));

        for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
                     cell_iter != rCellPopulation.End();
                     ++cell_iter)
        {
           boost::shared_ptr<AbstractCellProliferativeType> p_cell_type = cell_iter->GetCellProliferativeType();

           if(p_cell_type->IsType<StemCellProliferativeType>())
           {
                unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
                Node<DIM>* p_node = rCellPopulation.GetNode(node_index);

                double height = p_node->rGetLocation()[DIM-1];

                c_vector<double,DIM> force_direction = zero_vector<double>(DIM);
                force_direction[DIM-1]=-1.0;

                c_vector<double, DIM> force_contribution =mStemCellForceMagnitudeParameter*force_direction*(1.0+height*height);
                p_node->AddAppliedForceContribution(force_contribution);
            }
        }
    }
}

template<unsigned DIM>
void StemCellRetainerForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<StemCellForceMagnitudeParameter>" << mStemCellForceMagnitudeParameter << "</StemCellForceMagnitudeParameter> \n";

    // Call direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class StemCellRetainerForce<1>;
template class StemCellRetainerForce<2>;
template class StemCellRetainerForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(StemCellRetainerForce)
