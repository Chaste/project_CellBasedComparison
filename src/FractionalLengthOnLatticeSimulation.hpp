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

#ifndef FRACTIONALLENGTHONLATTICESIMULATION_HPP_
#define FRACTIONALLENGTHONLATTICESIMULATION_HPP_

#include <map>
#include "ChasteSerialization.hpp"
#include "OnLatticeSimulation.hpp"
#include "PetscTools.hpp"

/**
 * Subclass of OnLatticeSimulation to calculate and output the fractional length (the perimeter around labelled cells)
 */
template<unsigned DIM>
class FractionalLengthOnLatticeSimulation : public OnLatticeSimulation<DIM>
{
private :

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<OnLatticeSimulation<DIM> >(*this);
    }

    /**
     * Overridden SetupSolve() method. Calls UpdateCellData().
     */
    void SetupSolve();

    /**
     * Overridden UpdateAtEndOfTimeStep() method.
     */
    void UpdateAtEndOfTimeStep();

    /**
     * Overridden UpdateAtEndOfSolve() method.
     */
    void UpdateAtEndOfSolve();

    /**
     * Helper Method to Compute the fractional length and output to file
     */
    void CalculateFractionalLength();

public:

    /** The file to which the fractional lengths are written to. */
    out_stream mFractionalLengthsResultsFile;

    /**
     * Default constructor.
     *
     * @param rCellPopulation A cell population facade class (contains a mesh and cells)
     * @param deleteCellPopulationInDestructor Whether to delete the cell population on destruction to
     *     free up memory (defaults to false)
     * @param initialiseCells whether to initialise cells (defaults to true, set to false when loading from an archive)
     */
     FractionalLengthOnLatticeSimulation(AbstractCellPopulation<DIM>& rCellPopulation,
                                   bool deleteCellPopulationInDestructor=false,
                                   bool initialiseCells=true);

     /**
      * Destructor.
      */
    ~FractionalLengthOnLatticeSimulation();
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(FractionalLengthOnLatticeSimulation)

namespace boost
{
namespace serialization
{
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const FractionalLengthOnLatticeSimulation<DIM> * t, const BOOST_PFTO unsigned int file_version)
{
    const AbstractCellPopulation<DIM> * p_cell_population = &(t->rGetCellPopulation());
    ar & p_cell_population;
}

template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, FractionalLengthOnLatticeSimulation<DIM> * t, const unsigned int file_version)
{
    AbstractCellPopulation<DIM>* p_cell_population;
    ar >> p_cell_population;

    ::new(t)FractionalLengthOnLatticeSimulation<DIM>(*p_cell_population, true, false);
}
}
} // namespace ...

#endif /*FRACTIONALLENGTHONLATTICESIMULATION_HPP_*/
