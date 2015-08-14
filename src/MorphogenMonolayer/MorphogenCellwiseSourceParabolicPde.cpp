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

#include "MorphogenCellwiseSourceParabolicPde.hpp"

#include "AbstractCentreBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "CaBasedCellPopulation.hpp"
#include "ApoptoticCellProperty.hpp"
#include "Exception.hpp"
#include "CellLabel.hpp"
#include "Debug.hpp"

template<unsigned DIM>
MorphogenCellwiseSourceParabolicPde<DIM>::MorphogenCellwiseSourceParabolicPde(AbstractCellPopulation<DIM,DIM>& rCellPopulation,
        double duDtCoefficient,
        double diffusionCoefficient,
        double uptakeCoefficient,
        double sourceWidth)
		: CellwiseSourceParabolicPde<DIM>(rCellPopulation,duDtCoefficient, diffusionCoefficient, uptakeCoefficient ),
		  mSourceWidth(sourceWidth)
{
}

template<unsigned DIM>
const AbstractCellPopulation<DIM,DIM>& MorphogenCellwiseSourceParabolicPde<DIM>::rGetCellPopulation() const
{
    return this->mrCellPopulation;
}

template<unsigned DIM>
double MorphogenCellwiseSourceParabolicPde<DIM>::ComputeDuDtCoefficientFunction(const ChastePoint<DIM>& )
{
    return this->mDuDtCoefficient;
}


template<unsigned DIM>
double MorphogenCellwiseSourceParabolicPde<DIM>::ComputeSourceTerm(const ChastePoint<DIM>& rX, double u, Element<DIM,DIM>* pElement)
{
    NEVER_REACHED;
    return 0.0;
}

template<unsigned DIM>
c_matrix<double,DIM,DIM> MorphogenCellwiseSourceParabolicPde<DIM>::ComputeDiffusionTerm(const ChastePoint<DIM>& rX, Element<DIM,DIM>* pElement)
{
    return this->mDiffusionCoefficient*identity_matrix<double>(DIM);
}

template<unsigned DIM>
double MorphogenCellwiseSourceParabolicPde<DIM>::ComputeSourceTermAtNode(const Node<DIM>& rNode, double u)
{
    double coefficient = 0.0;

    double x = rNode.rGetLocation()[0];
    //double y = rNode.rGetLocation()[1];

    //if (x*x+y*y <4)

    if (x>-0.5*mSourceWidth && x<0.5*mSourceWidth)
    {
    	coefficient = this->mUptakeCoefficient;
    }

    double decay_coeficient = 0.01;

    return coefficient - u*decay_coeficient;
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class MorphogenCellwiseSourceParabolicPde<1>;
template class MorphogenCellwiseSourceParabolicPde<2>;
template class MorphogenCellwiseSourceParabolicPde<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(MorphogenCellwiseSourceParabolicPde)
