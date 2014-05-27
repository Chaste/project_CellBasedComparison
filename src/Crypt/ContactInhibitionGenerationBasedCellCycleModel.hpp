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

#ifndef ContactInhibitionGenerationBasedCellCycleModel_HPP_
#define ContactInhibitionGenerationBasedCellCycleModel_HPP_

#include "AbstractSimpleGenerationBasedCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"

/**
 * A stochastic cell-cycle model employed by Meineke et al (2001) in their off-lattice
 * model of the intestinal crypt (doi:10.1046/j.0960-7722.2001.00216.x). WIth Added Contact Inhibition
 */
class ContactInhibitionGenerationBasedCellCycleModel : public AbstractSimpleGenerationBasedCellCycleModel
{
    friend class TestSimpleCellCycleModels;

private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the cell-cycle model and random number generator, never used directly - boost uses this.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractSimpleGenerationBasedCellCycleModel>(*this);
        archive & mQuiescentVolumeFraction;
        archive & mEquilibriumVolume;
        archive & mCurrentQuiescentDuration;
        archive & mCurrentQuiescentOnsetTime;

        // Make sure the RandomNumberGenerator singleton gets saved too
        SerializableSingleton<RandomNumberGenerator>* p_wrapper = RandomNumberGenerator::Instance()->GetSerializationWrapper();
        archive & p_wrapper;
    }

    /**
	* The fraction of the cells' equilibrium volume in G1 phase below which these cells are quiescent.
	*/
	double mQuiescentVolumeFraction;

	/**
	* The cell equilibrium volume while in G1 phase.
	*/
	double mEquilibriumVolume;

	/**
	* The time when the current period of quiescence began.
	*/
	double mCurrentQuiescentOnsetTime;

	/**
	* How long the current period of quiescence has lasted.
	* Has units of hours.
	*/
	double mCurrentQuiescentDuration;

    /**
     * Stochastically set the G1 duration.  Called on cell creation at
     * the start of a simulation, and for both parent and daughter
     * cells at cell division.
     */
    void SetG1Duration();

public:

    /**
     * Constructor - just a default, mBirthTime is now set in the AbstractCellCycleModel class.
     * mG1Duration is set very high, it is set for the individual cells when InitialiseDaughterCell is called
     */
    ContactInhibitionGenerationBasedCellCycleModel();

    /**
     * Overridden builder method to create new copies of
     * this cell-cycle model.
     */
    AbstractCellCycleModel* CreateCellCycleModel();

    /**
     * Overridden UpdateCellCyclePhase() method to include Contact inhibition.
     */
    virtual void UpdateCellCyclePhase();

    /**
     * @param quiescentVolumeFraction
     */
    void SetQuiescentVolumeFraction(double quiescentVolumeFraction);

    /**
     * @return mQuiescentVolumeFraction
     */
    double GetQuiescentVolumeFraction();

    /**
     * @param equilibriumVolume
     */
    void SetEquilibriumVolume(double equilibriumVolume);

    /**
     * @return mEquilibriumVolume
     */
    double GetEquilibriumVolume();

    /**
    * Set method for mCurrentQuiescentDuration.
    *
    * @param currentQuiescentDuration the new value of mCurrentQuiescentDuration
    */
    void SetCurrentQuiescentDuration(double currentQuiescentDuration);

    /**
     * @return mCurrentQuiescentDuration
     */
    double GetCurrentQuiescentDuration();

    /**
    * Set method for mCurrentQuiescentOnsetTime.
    *
    * @param currentQuiescentOnsetTime the new value of mCurrentQuiescentOnsetTime
    */
    void SetCurrentQuiescentOnsetTime(double currentQuiescentOnsetTime);

    /**
     * @return mCurrentQuiescentOnsetTime
     */
    double GetCurrentQuiescentOnsetTime();

    /**
     * Outputs cell cycle model parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputCellCycleModelParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(ContactInhibitionGenerationBasedCellCycleModel)

#endif /*ContactInhibitionGenerationBasedCellCycleModel_HPP_*/
