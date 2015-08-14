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

#ifndef CELLWISESOURCEMORPHOGENPDE_HPP_
#define CELLWISESOURCEMORPHOGENPDE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>


#include "CellwiseSourceParabolicPde.hpp"

/**
 *  A PDE which has a source at each labelled non-apoptotic cell.
 *
 *  Modified from CellwiseSourceParabolicPde so only labelled cells have a source
 *
 */
template<unsigned DIM>
class MorphogenCellwiseSourceParabolicPde : public CellwiseSourceParabolicPde<DIM>
{

private:
    /*
     * Stores how side the central region of source cells is.
     */
    double mSourceWidth;

    /** Needed for serialization.*/
    friend class boost::serialization::access;
    /**
     * Serialize the PDE and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
       archive & boost::serialization::base_object<CellwiseSourceParabolicPde<DIM> >(*this);
       archive & mSourceWidth;
    }

public:

    /**
	 * Constructor.
	 *
	 * @param rCellPopulation reference to the cell population
	 * @param duDtCoefficient rate of reaction (defaults to 1.0)
	 * @param diffusionCoefficient rate of diffusion (defaults to 1.0)
	 * @param uptakeCoefficient the coefficient of consumption of nutrient by cells (defaults to 0.0)
	 * @param sourceWidth the width of the source of morphgen in the centte of the mesh (defaults to 1.0)
	 */
    MorphogenCellwiseSourceParabolicPde(AbstractCellPopulation<DIM, DIM>& rCellPopulation,
                                   double duDtCoefficient = 1.0,
                                   double diffusionCoefficient = 1.0,
                                   double uptakeCoefficient = 0.0,
                                   double sourceWidth = 2.0);

    /**
         * @return const reference to the cell population (used in archiving).
         */
        const AbstractCellPopulation<DIM>& rGetCellPopulation() const;


    /**
	 * Overridden ComputeDuDtCoefficientFunction() method
	 *
	 * @return the function c(x) in "c(x) du/dt = Grad.(DiffusionTerm(x)*Grad(u))+LinearSourceTerm(x)+NonlinearSourceTerm(x, u)"
	 *
	 * @param rX the point in space at which the function c is computed
	 */
	virtual double ComputeDuDtCoefficientFunction(const ChastePoint<DIM>& rX);

	/**
	 * Overridden ComputeSourceTerm() method.
	 *
	 * @return computed source term.
	 *
	 * @param rX the point in space at which the nonlinear source term is computed
	 * @param u the value of the dependent variable at the point
	 */
	virtual double ComputeSourceTerm(const ChastePoint<DIM>& rX,
									 double u,
									 Element<DIM,DIM>* pElement=NULL);


	/**
	 * Overridden ComputeSourceTermAtNode() method.
	 *
	 * @return computed source term at a node.
	 *
	 * @param rNode the node at which the nonlinear source term is computed
	 * @param u the value of the dependent variable at the node
	 */
	virtual double ComputeSourceTermAtNode(const Node<DIM>& rNode, double u);

	/**
	 * Overridden ComputeDiffusionTerm() method.
	 *
	 * @param rX The point in space at which the diffusion term is computed
	 * @param pElement The mesh element that x is contained in (optional).
	 *
	 * @return a matrix.
	 */
	virtual c_matrix<double,DIM,DIM> ComputeDiffusionTerm(const ChastePoint<DIM>& rX, Element<DIM,DIM>* pElement=NULL);

};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(MorphogenCellwiseSourceParabolicPde)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a MorphogenCellwiseSourceParabolicPde.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const MorphogenCellwiseSourceParabolicPde<DIM>* t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<DIM, DIM>* p_cell_population = &(t->rGetCellPopulation());
    ar & p_cell_population;
}

/**
 * De-serialize constructor parameters and initialise a MorphogenCellwiseSourceParabolicPde.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, MorphogenCellwiseSourceParabolicPde<DIM>* t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<DIM, DIM>* p_cell_population;
    ar >> p_cell_population;

    // Invoke inplace constructor to initialise instance
    ::new(t)MorphogenCellwiseSourceParabolicPde<DIM>(*p_cell_population);
}
}
} // namespace ...

#endif /*CELLWISESOURCEMORPHOGENPDE_HPP_*/
