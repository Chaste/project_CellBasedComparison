/*

Copyright (c) 2005-2014, University of Oxford.
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

#include "StemCellRetainerPottsUpdateRule.hpp"

template<unsigned DIM>
StemCellRetainerPottsUpdateRule<DIM>::StemCellRetainerPottsUpdateRule()
    : AbstractPottsUpdateRule<DIM>(),
      mStemCellRestraintParameter(0.1)
{
}

template<unsigned DIM>
StemCellRetainerPottsUpdateRule<DIM>::~StemCellRetainerPottsUpdateRule()
{
}

template<unsigned DIM>
double StemCellRetainerPottsUpdateRule<DIM>::EvaluateHamiltonianContribution(unsigned currentNodeIndex,
                                                                        unsigned targetNodeIndex,
                                                                        PottsBasedCellPopulation<DIM>& rCellPopulation)
{
    double delta_H = 0.0;

    //See if the target node is part of a stem cell
    std::set<unsigned> new_location_containing_elements = rCellPopulation.GetNode(targetNodeIndex)->rGetContainingElementIndices();
    bool target_node_contained = !new_location_containing_elements.empty();

    if (target_node_contained)
    {
    	unsigned target_element_index = (*new_location_containing_elements.begin());

    	CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(target_element_index);

		if (p_cell->GetCellProliferativeType()->IsType<StemCellProliferativeType>())
		{
			// Note that we define these vectors before setting them as otherwise the profiling build will break (see #2367)
			c_vector<double, DIM> current_location;
			current_location = rCellPopulation.GetNode(currentNodeIndex)->rGetLocation();
			c_vector<double, DIM> target_location;
			target_location = rCellPopulation.GetNode(targetNodeIndex)->rGetLocation();

			if (target_location[DIM-1] > current_location[DIM-1])
			{
				delta_H += mStemCellRestraintParameter;
			}
			else if (target_location[DIM-1] < current_location[DIM-1])
			{
				delta_H -= mStemCellRestraintParameter;
			}
		}
    }
    return delta_H;
}

template<unsigned DIM>
void StemCellRetainerPottsUpdateRule<DIM>::SetStemCellRestraintParameter(double stemCellRestraintParameter)
{
	mStemCellRestraintParameter = stemCellRestraintParameter;
}

template<unsigned DIM>
double StemCellRetainerPottsUpdateRule<DIM>::GetStemCellRestraintParameter()
{
    return mStemCellRestraintParameter;
}

template<unsigned DIM>
void StemCellRetainerPottsUpdateRule<DIM>::OutputUpdateRuleParameters(out_stream& rParamsFile)
{
	  *rParamsFile << "\t\t\t<StemCellRestraintParameter>" << mStemCellRestraintParameter << "</StemCellRestraintParameter> \n";

    // Call method on direct parent class
    AbstractPottsUpdateRule<DIM>::OutputUpdateRuleParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class StemCellRetainerPottsUpdateRule<1>;
template class StemCellRetainerPottsUpdateRule<2>;
template class StemCellRetainerPottsUpdateRule<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(StemCellRetainerPottsUpdateRule)
