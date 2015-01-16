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

#include "PeriodicNodeBasedBoundaryCondition.hpp"
#include "WntConcentration.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"
#include "RandomNumberGenerator.hpp"
#include "Debug.hpp"

PeriodicNodeBasedBoundaryCondition::PeriodicNodeBasedBoundaryCondition(AbstractCellPopulation<2>* pCellPopulation, double width)
    : AbstractCellPopulationBoundaryCondition<2>(pCellPopulation),
      mWidth(width)
{
}

void PeriodicNodeBasedBoundaryCondition::ImposeBoundaryCondition(const std::map<Node<2>*, c_vector<double, 2> >& rOldLocations)
{
    assert (dynamic_cast<AbstractCentreBasedCellPopulation<2>*>(this->mpCellPopulation));



    // Iterate over all nodes to update their positions. Using iterator to ignore deleted nodes.
    for (AbstractMesh<2,2>::NodeIterator node_iter = this->mpCellPopulation->rGetMesh().GetNodeIteratorBegin();
             node_iter != this->mpCellPopulation->rGetMesh().GetNodeIteratorEnd();
             ++node_iter)
    {
        // Any node that has moved below the bottom of the crypt must be moved back up
        double x = node_iter->rGetLocation()[0];
        if (x < 0.0)
        {
            node_iter->rGetModifiableLocation()[0] = x + mWidth;
        }
        else if (x > mWidth)
        {
            node_iter->rGetModifiableLocation()[0] = x - mWidth;
        }

        assert(node_iter->rGetLocation()[0] >= 0.0);
        assert(node_iter->rGetLocation()[0] <= mWidth);
    }
}

bool PeriodicNodeBasedBoundaryCondition::VerifyBoundaryCondition()
{
    bool boundary_condition_satisfied = true;

    return boundary_condition_satisfied;
}


void PeriodicNodeBasedBoundaryCondition::SetWidth(bool width)
{
    mWidth = width;
}

bool PeriodicNodeBasedBoundaryCondition::GetWidth()
{
    return mWidth;
}

void PeriodicNodeBasedBoundaryCondition::OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t<Width>" << mWidth << "</Width>\n";

    // Call method on direct parent class
    AbstractCellPopulationBoundaryCondition<2>::OutputCellPopulationBoundaryConditionParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(PeriodicNodeBasedBoundaryCondition)
