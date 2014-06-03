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

#include "DeltaNotchWriter.hpp"
#include "AbstractCellPopulation.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
DeltaNotchWriter<ELEMENT_DIM, SPACE_DIM>::DeltaNotchWriter()
    : AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("celldeltanotch.dat")
{
    this->mVtkCellDataName = "delta";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double DeltaNotchWriter<ELEMENT_DIM, SPACE_DIM>::GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    double delta= pCell->GetCellData()->GetItem("delta");
    return delta;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void DeltaNotchWriter<ELEMENT_DIM, SPACE_DIM>::VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    unsigned location_index = pCellPopulation->GetLocationIndexUsingCell(pCell);
    unsigned cell_id = pCell->GetCellId();
    c_vector<double, SPACE_DIM> centre_location = pCellPopulation->GetLocationOfCellCentre(pCell);
    double delta= pCell->GetCellData()->GetItem("delta");
    double notch= pCell->GetCellData()->GetItem("notch");
    double mean_delta= pCell->GetCellData()->GetItem("mean delta");

	*this->mpOutStream << location_index << " " << cell_id << " ";
	for (unsigned i=0; i<SPACE_DIM; i++)
	{
		*this->mpOutStream << centre_location[i] << " ";
	}

	*this->mpOutStream << delta << " " << notch << " " << mean_delta << " ";
}

// Explicit instantiation
template class DeltaNotchWriter<1,1>;
template class DeltaNotchWriter<1,2>;
template class DeltaNotchWriter<2,2>;
template class DeltaNotchWriter<1,3>;
template class DeltaNotchWriter<2,3>;
template class DeltaNotchWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(DeltaNotchWriter)
