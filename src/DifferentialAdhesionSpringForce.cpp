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

#include "DifferentialAdhesionSpringForce.hpp"
#include "CellLabel.hpp"
#include "Debug.hpp"

template<unsigned DIM>
DifferentialAdhesionSpringForce<DIM>::DifferentialAdhesionSpringForce()
   : GeneralisedLinearSpringForce<DIM>(),
     mDifferentialAdhesionFactor(0.5)
{
}

template<unsigned DIM>
double DifferentialAdhesionSpringForce<DIM>::VariableSpringConstantMultiplicationFactor(
    unsigned nodeAGlobalIndex,
    unsigned nodeBGlobalIndex,
    AbstractCellPopulation<DIM>& rCellPopulation,
    bool isCloserThanRestLength)
{
	CellPtr p_cell_A = rCellPopulation.GetCellUsingLocationIndex(nodeAGlobalIndex);
    CellPtr p_cell_B = rCellPopulation.GetCellUsingLocationIndex(nodeBGlobalIndex);

	// Less attraction between different types.
	if (( p_cell_A->template HasCellProperty<CellLabel>()  &&
		  p_cell_B->template HasCellProperty<CellLabel>() ) ||
		(!(p_cell_A->template HasCellProperty<CellLabel>())  &&
		 !(p_cell_B->template HasCellProperty<CellLabel>()) ) )
	{
		return 1.0;
	}
	else if ((!(p_cell_A->template HasCellProperty<CellLabel>())  &&
			   (p_cell_B->template HasCellProperty<CellLabel>()) ) ||
			 ( (p_cell_A->template HasCellProperty<CellLabel>())  &&
			  !(p_cell_B->template HasCellProperty<CellLabel>()) ) )
	{
		if (isCloserThanRestLength)
		{
			return mDifferentialAdhesionFactor; //Repulsion
		}
		else
		{
			return mDifferentialAdhesionFactor; // Attraction
		}
	}
	else
	{
		NEVER_REACHED;
		return 0.0;
	}
}

template<unsigned DIM>
c_vector<double, DIM> DifferentialAdhesionSpringForce<DIM>::CalculateForceBetweenNodes(unsigned nodeAGlobalIndex,
                                                                                    unsigned nodeBGlobalIndex,
                                                                                    AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // We should only ever calculate the force between two distinct nodes
    assert(nodeAGlobalIndex != nodeBGlobalIndex);

    Node<DIM>* p_node_a = rCellPopulation.GetNode(nodeAGlobalIndex);
    Node<DIM>* p_node_b = rCellPopulation.GetNode(nodeBGlobalIndex);

    // Get the node locations
    c_vector<double, DIM> node_a_location = p_node_a->rGetLocation();
    c_vector<double, DIM> node_b_location = p_node_b->rGetLocation();

    // Get the node radii for a NodeBasedCellPopulation
    double node_a_radius=0.0;
    double node_b_radius=0.0;

    if (dynamic_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation))
    {
        node_a_radius = p_node_a->GetRadius();
        node_b_radius = p_node_b->GetRadius();
    }

    // Get the unit vector parallel to the line joining the two nodes
    c_vector<double, DIM> unit_difference;
    /*
     * We use the mesh method GetVectorFromAtoB() to compute the direction of the
     * unit vector along the line joining the two nodes, rather than simply subtract
     * their positions, because this method can be overloaded (e.g. to enforce a
     * periodic boundary in Cylindrical2dMesh).
     */
    unit_difference = rCellPopulation.rGetMesh().GetVectorFromAtoB(node_a_location, node_b_location);

    // Calculate the distance between the two nodes
    double distance_between_nodes = norm_2(unit_difference);
    assert(distance_between_nodes > 0);
    assert(!std::isnan(distance_between_nodes));

    unit_difference /= distance_between_nodes;

    /*
     * If mUseCutOffLength has been set, then there is zero force between
     * two nodes located a distance apart greater than mMechanicsCutOffLength in AbstractTwoBodyInteractionForce.
     */
    if (this->mUseCutOffLength)
    {
        if (distance_between_nodes >= this->GetCutOffLength())
        {
            return zero_vector<double>(DIM); // c_vector<double,DIM>() is not guaranteed to be fresh memory
        }
    }

    /*
     * Use a default rest length of 1.0 for all springs
     */
    double rest_length = 1.0;


    double overlap = distance_between_nodes - rest_length;
    bool is_closer_than_rest_length = (overlap <= 0);
    double multiplication_factor = VariableSpringConstantMultiplicationFactor(nodeAGlobalIndex, nodeBGlobalIndex, rCellPopulation, is_closer_than_rest_length);
    double spring_stiffness = this->mMeinekeSpringStiffness;

    {
        // A reasonably stable simple force law
        if (is_closer_than_rest_length) //overlap is negative
        {
            //log(x+1) is undefined for x<=-1
            assert(overlap > -	rest_length);
            c_vector<double, DIM> temp = multiplication_factor*spring_stiffness * unit_difference * rest_length* log(1.0 + overlap/rest_length);
            return temp;
        }
        else
        {
            double alpha = 4;
            c_vector<double, DIM> temp = multiplication_factor*spring_stiffness * unit_difference * overlap * exp(-alpha * overlap/rest_length);
            return temp;
        }
    }
}


template<unsigned DIM>
double DifferentialAdhesionSpringForce<DIM>::GetDifferentialAdhesionFactor()
{
    return mDifferentialAdhesionFactor;
}

template<unsigned DIM>
void DifferentialAdhesionSpringForce<DIM>::SetDifferentialAdhesionFactor(double differentialAdhesionFactor)
{
    assert(differentialAdhesionFactor > 0.0);
    mDifferentialAdhesionFactor = differentialAdhesionFactor;
}


template<unsigned DIM>
void DifferentialAdhesionSpringForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
	*rParamsFile << "\t\t\t<DifferentialAdhesionFactor>" << mDifferentialAdhesionFactor << "</DifferentialAdhesionFactor>\n";

	// Call direct parent class
    GeneralisedLinearSpringForce<DIM>::OutputForceParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class DifferentialAdhesionSpringForce<1>;
template class DifferentialAdhesionSpringForce<2>;
template class DifferentialAdhesionSpringForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(DifferentialAdhesionSpringForce)
