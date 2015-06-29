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

#include "RandomMotionForce.hpp"

template<unsigned DIM>
RandomMotionForce<DIM>::RandomMotionForce()
    : AbstractForce<DIM>(),
	  mMovementParameter(0.01)
{
}

template<unsigned DIM>
RandomMotionForce<DIM>::~RandomMotionForce()
{
}

template<unsigned DIM>
void RandomMotionForce<DIM>::SetMovementParameter(double movementParameter)
{
    assert(movementParameter > 0.0);
    mMovementParameter = movementParameter;
}

template<unsigned DIM>
double RandomMotionForce<DIM>::GetMovementParameter()
{
    return mMovementParameter;
}

template<unsigned DIM>
void RandomMotionForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
    double dt = SimulationTime::Instance()->GetTimeStep();

    // Loop over the cells
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        // Get the node index associated with this cell
        unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
        Node<DIM>* p_node = rCellPopulation.GetNode(node_index);

        c_vector<double, DIM> force_contribution;
        for (unsigned i=0; i<DIM; i++)
        {
            /*
             * The force on this cell is scaled with the timestep such that when it is
             * used in the discretised equation of motion for the cell, we obtain the
             * correct formula
             *
             * x_new = x_old + sqrt(2*D*dt)*W
             *
             * where W is a standard normal random variable.
             */
            double xi = RandomNumberGenerator::Instance()->StandardNormalRandomDeviate();

            force_contribution[i] = (sqrt(2.0*mMovementParameter*dt)/dt)*xi;
        }
        p_node->AddAppliedForceContribution(force_contribution);
    }
}

template<unsigned DIM>
void RandomMotionForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<MovementParameter>" << mMovementParameter << "</MovementParameter> \n";

    // Call direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class RandomMotionForce<1>;
template class RandomMotionForce<2>;
template class RandomMotionForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(RandomMotionForce)
