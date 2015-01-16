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

#include "ContactInhibitionGenerationBasedCellCycleModel.hpp"
#include "Exception.hpp"
#include "CellLabel.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "PanethCellProliferativeType.hpp"

ContactInhibitionGenerationBasedCellCycleModel::ContactInhibitionGenerationBasedCellCycleModel()
    : AbstractSimpleGenerationBasedCellCycleModel(),
      mQuiescentVolumeFraction(DOUBLE_UNSET),
     mEquilibriumVolume(DOUBLE_UNSET),
     mCurrentQuiescentOnsetTime(SimulationTime::Instance()->GetTime()),
     mCurrentQuiescentDuration(0.0)
{
}

AbstractCellCycleModel* ContactInhibitionGenerationBasedCellCycleModel::CreateCellCycleModel()
{
    // Create a new cell-cycle model
    ContactInhibitionGenerationBasedCellCycleModel* p_model = new ContactInhibitionGenerationBasedCellCycleModel();

    /*
     * Set each member variable of the new cell-cycle model that inherits
     * its value from the parent.
     *
     * Note 1: some of the new cell-cycle model's member variables (namely
     * mBirthTime, mCurrentCellCyclePhase, mReadyToDivide) will already have been
     * correctly initialized in its constructor.
     *
     * Note 2: one or more of the new cell-cycle model's member variables
     * may be set/overwritten as soon as InitialiseDaughterCell() is called on
     * the new cell-cycle model.
     *
     * Note 3: the member variable mDimension remains unset, since this cell-cycle
     * model does not need to know the spatial dimension, so if we were to call
     * SetDimension() on the new cell-cycle model an exception would be triggered;
     * hence we do not set this member variable.
     */
    p_model->SetBirthTime(mBirthTime);
    p_model->SetMinimumGapDuration(mMinimumGapDuration);
    p_model->SetStemCellG1Duration(mStemCellG1Duration);
    p_model->SetTransitCellG1Duration(mTransitCellG1Duration);
    p_model->SetSDuration(mSDuration);
    p_model->SetG2Duration(mG2Duration);
    p_model->SetMDuration(mMDuration);
    p_model->SetGeneration(mGeneration);
    p_model->SetMaxTransitGenerations(mMaxTransitGenerations);
    p_model->SetQuiescentVolumeFraction(mQuiescentVolumeFraction);
    p_model->SetEquilibriumVolume(mEquilibriumVolume);
    p_model->SetCurrentQuiescentOnsetTime(mCurrentQuiescentOnsetTime);
    p_model->SetCurrentQuiescentDuration(mCurrentQuiescentDuration);

    return p_model;
}

void ContactInhibitionGenerationBasedCellCycleModel::UpdateCellCyclePhase()
{
 /// Copied from   ContactInhibitionBasedCellCycleModel ////
    // Get cell volume
    double cell_volume = mpCell->GetCellData()->GetItem("volume");

    // Removes the cell label
    mpCell->RemoveCellProperty<CellLabel>();

    if (mCurrentCellCyclePhase == G_ONE_PHASE)
    {
        // Update G1 duration based on cell volume
        double dt = SimulationTime::Instance()->GetTimeStep();
        double quiescent_volume = mEquilibriumVolume * mQuiescentVolumeFraction;

        if (cell_volume < quiescent_volume)
        {
            // Update the duration of the current period of contact inhibition.
            mCurrentQuiescentDuration = SimulationTime::Instance()->GetTime() - mCurrentQuiescentOnsetTime;
            mG1Duration += dt;
            mpCell->AddCellProperty(CellPropertyRegistry::Instance()->Get<CellLabel>());
        }
         else
        {
            // Reset the cell's quiescent duration and update the time at which the onset of quiescent occurs
            mCurrentQuiescentDuration = 0.0;
            mCurrentQuiescentOnsetTime = SimulationTime::Instance()->GetTime();
        }
    }
////////////////////////////////////////////////////////

    double time_since_birth = GetAge();
    assert(time_since_birth >= 0);

    if (mpCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>()
            || mpCell->GetCellProliferativeType()->IsType<PanethCellProliferativeType>())
    {
        mCurrentCellCyclePhase = G_ZERO_PHASE;
    }
    else if (time_since_birth < GetMDuration())
    {
        mCurrentCellCyclePhase = M_PHASE;
    }
    else if (time_since_birth < GetMDuration() + mG1Duration)
    {
        mCurrentCellCyclePhase = G_ONE_PHASE;
    }
    else if (time_since_birth < GetMDuration() + mG1Duration + GetSDuration())
    {
        mCurrentCellCyclePhase = S_PHASE;
    }
    else if (time_since_birth < GetMDuration() + mG1Duration + GetSDuration() + GetG2Duration())
    {
        mCurrentCellCyclePhase = G_TWO_PHASE;
    }
}

void ContactInhibitionGenerationBasedCellCycleModel::SetG1Duration()
{
    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();

    if (mpCell->GetCellProliferativeType()->IsType<StemCellProliferativeType>())
    {
        mG1Duration = GetStemCellG1Duration() + 4*p_gen->ranf(); // U[14,18] for default parameters (mStemCellG1Duration) according to Meineke
    }
    else if (mpCell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>())
    {
        mG1Duration = GetTransitCellG1Duration() + 2*p_gen->ranf(); // U[4,6] for default parameters (mTransitG1CellDuration) according to Meineke
    }
    else if (mpCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>()
            || mpCell->GetCellProliferativeType()->IsType<PanethCellProliferativeType>())
    {
       mG1Duration = DBL_MAX;
    }
    else
    {
       NEVER_REACHED;
    }


}

void ContactInhibitionGenerationBasedCellCycleModel::SetQuiescentVolumeFraction(double quiescentVolumeFraction)
{
    mQuiescentVolumeFraction = quiescentVolumeFraction;
}

double ContactInhibitionGenerationBasedCellCycleModel::GetQuiescentVolumeFraction()
{
    return mQuiescentVolumeFraction;
}

void ContactInhibitionGenerationBasedCellCycleModel::SetEquilibriumVolume(double equilibriumVolume)
{
    mEquilibriumVolume = equilibriumVolume;
}

double ContactInhibitionGenerationBasedCellCycleModel::GetEquilibriumVolume()
{
    return mEquilibriumVolume;
}

void ContactInhibitionGenerationBasedCellCycleModel::SetCurrentQuiescentDuration(double currentQuiescentDuration)
{
    mCurrentQuiescentDuration = currentQuiescentDuration;
}

double ContactInhibitionGenerationBasedCellCycleModel::GetCurrentQuiescentDuration()
{
    return mCurrentQuiescentDuration;
}

void ContactInhibitionGenerationBasedCellCycleModel::SetCurrentQuiescentOnsetTime(double currentQuiescentOnsetTime)
{
    mCurrentQuiescentOnsetTime = currentQuiescentOnsetTime;
}

double ContactInhibitionGenerationBasedCellCycleModel::GetCurrentQuiescentOnsetTime()
{
    return mCurrentQuiescentOnsetTime;
}

void ContactInhibitionGenerationBasedCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<QuiescentVolumeFraction>" << mQuiescentVolumeFraction << "</QuiescentVolumeFraction>\n";
    *rParamsFile << "\t\t\t<EquilibriumVolume>" << mEquilibriumVolume << "</EquilibriumVolume>\n";

    // Call method on direct parent class
    AbstractSimpleGenerationBasedCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(ContactInhibitionGenerationBasedCellCycleModel)
