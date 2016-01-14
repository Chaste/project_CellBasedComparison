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

#ifndef MORPHOGENDEPENDENTCELLCYCLEMODEL_HPP_
#define MORPHOGENDEPENDENTCELLCYCLEMODEL_HPP_

#include "AbstractCellCycleModel.hpp"

/**
 * TODO
 */

class MorphogenDependentCellCycleModel : public AbstractCellCycleModel
{
private:

    friend class boost::serialization::access;

    /**
     * Boost Serialization method for archiving/checkpointing
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellCycleModel>(*this);
        archive & mGrowthRate;
        archive & mCurrentMass;
        archive & mTargetMass;
        archive & mMorphogenInfluence;
    }

    /**
     * Growth rate of the cell, randomly assigned on birth
     */
    double mGrowthRate;

    /**
     * Stores the current mass of the cell, (Referred to as volume in the paper)
     */
    double mCurrentMass;

    // Parameters

    /**
     * Target mass of the cell (Referred to as volume in the paper)
     */
    double mTargetMass;

    /**
     * the level of influence the morphogen has on growth
     */
    double mMorphogenInfluence;

public:

    /**
     * Constructor.
     */
    MorphogenDependentCellCycleModel();

    /**
     * Overridden ReadyToDivideMethod
     *
     * @return whether the cell is ready to divide (enter M phase).
     *
     * Here we divide randomly basd on the size of the cell compared to the target size.
     *
     */
    virtual bool ReadyToDivide();

    /**
     * Overriden ResetForDivision method to halve the current cell size.
     *
     * Should only be called by the Cell Divide() method.
     */
    virtual void ResetForDivision();

    /**
     * Overridden UpdateCellCyclePhase() method.
     *
     * Note this is never called. But we need to include it.
     */
    void UpdateCellCyclePhase();

    /**
     * Overridden builder method to create new instances of
     * the cell-cycle model.
     *
     * @return new cell-cycle model
     *
     */
    AbstractCellCycleModel* CreateCellCycleModel();

    /**
     * @param thresholdMorphogen
     */
    void SetCurrentMass(double currentMass);

    /**
     * @return mCurrentMass
     */
    double GetCurrentMass();

    /**
     * @param targetMass
     */
    void SetTargetMass(double targetMass);

    /**
     * @return mTargetMass
     */
    double GetTargetMass();

    /**
     * @param morphogenInfluence
     */
    void SetMorphogenInfluence(double morphogenInfluence);

    /**
     * @return mMorphogenInfluence
     */
    double GetMorphogenInfluence();

    /**
     * Helper method to generate the mGrowthRate from a truncated normal distribution.
     */
    void GenerateGrowthRate();

    /**
     * Outputs cell cycle model parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputCellCycleModelParameters(out_stream& rParamsFile);
};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(MorphogenDependentCellCycleModel)

#endif // MORPHOGENDEPENDENTCELLCYCLEMODEL_HPP_
