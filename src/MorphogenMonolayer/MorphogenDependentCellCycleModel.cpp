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

#include "MorphogenDependentCellCycleModel.hpp"
#include "CellLabel.hpp"
#include "DefaultCellProliferativeType.hpp"
#include "RandomNumberGenerator.hpp"

MorphogenDependentCellCycleModel::MorphogenDependentCellCycleModel()
    : AbstractCellCycleModel(),
	  mGrowthRate(DOUBLE_UNSET),
	  mCurrentMass(DOUBLE_UNSET),
	  mTargetMass(1.0),
	  mMorphogenInfluence(DOUBLE_UNSET)
{

	double minimum_growth_rate = 0.01;
	double mean_growth_rate = 0.05;
	double sd = 0.01;
    RandomNumberGenerator* p_gen =RandomNumberGenerator::Instance();

	//Genrate a truncated normal for the Growth rate
	mGrowthRate = 0.0;
	while ( mGrowthRate < minimum_growth_rate)
	{
		mGrowthRate = p_gen->NormalRandomDeviate(mean_growth_rate,sd);

	}
}

bool MorphogenDependentCellCycleModel::ReadyToDivide()
{
	if  ((mCurrentMass == DOUBLE_UNSET) || (mMorphogenInfluence == DOUBLE_UNSET) )
    {
        EXCEPTION("The member variables mCurrentMass and mMorphogenInfluence have not yet been set.");
    }

	assert(mpCell != NULL);

	double morphogen_level = mpCell->GetCellData()->GetItem("morphogen");
	double dt = SimulationTime::Instance()->GetTimeStep();


	// Equation 16 from Smith et al 2011
	mCurrentMass *= 1.0 + dt*mGrowthRate*(1.0+mMorphogenInfluence*morphogen_level)*(1.0-mCurrentMass/mTargetMass);

	double division_rate = 0.1; // Per cell per hour


	double p_division = division_rate * dt *  (mCurrentMass/mTargetMass - 0.5);
	RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
	 if (p_gen->ranf()<p_division)
	 {
		 mReadyToDivide = true;
	 }

	 mpCell->GetCellData()->SetItem("p_division", p_division);
	 mpCell->GetCellData()->SetItem("GrowthRate", mGrowthRate);
	 mpCell->GetCellData()->SetItem("CurrentMass", mCurrentMass);

	return mReadyToDivide;
}

void MorphogenDependentCellCycleModel::ResetForDivision()
{
    AbstractCellCycleModel::ResetForDivision();
    mBirthTime = SimulationTime::Instance()->GetTime();
    // Halve the mass as half goes to each daughter cell
    mCurrentMass *= 0.5;
}

void MorphogenDependentCellCycleModel::UpdateCellCyclePhase()
{
	NEVER_REACHED;
}

AbstractCellCycleModel* MorphogenDependentCellCycleModel::CreateCellCycleModel()
{
    // Create a new cell-cycle model
    MorphogenDependentCellCycleModel* p_model = new MorphogenDependentCellCycleModel();

    /*
     * Set each member variable of the new cell-cycle model that inherits
     * its value from the parent.
     *
     * Note 1: some of the new cell-cycle model's member variables (namely
     * mBirthTime, mCurrentCellCyclePhase, mReadyToDivide, mTimeSpentInG1Phase,
     * mCurrentHypoxicDuration, mCurrentHypoxiaOnsetTime) will already have been
     * correctly initialized in its constructor.
     *
     * Note 2: one or more of the new cell-cycle model's member variables
     * may be set/overwritten as soon as InitialiseDaughterCell() is called on
     * the new cell-cycle model.
     */
    p_model->SetBirthTime(mBirthTime);
    p_model->SetDimension(mDimension);
    p_model->SetMinimumGapDuration(mMinimumGapDuration);
    p_model->SetStemCellG1Duration(mStemCellG1Duration);
    p_model->SetTransitCellG1Duration(mTransitCellG1Duration);
    p_model->SetSDuration(mSDuration);
    p_model->SetG2Duration(mG2Duration);
    p_model->SetMDuration(mMDuration);
    p_model->SetCurrentMass(mCurrentMass);
    p_model->SetTargetMass(mTargetMass);
    p_model->SetMorphogenInfluence(mMorphogenInfluence);

    return p_model;
}


void MorphogenDependentCellCycleModel::SetCurrentMass(double currentMass)
{
	mCurrentMass = currentMass;
}

double MorphogenDependentCellCycleModel::GetCurrentMass()
{
    return mCurrentMass;
}

void MorphogenDependentCellCycleModel::SetTargetMass(double targetMass)
{
	mTargetMass = targetMass;
}

double MorphogenDependentCellCycleModel::GetTargetMass()
{
    return mTargetMass;
}

void MorphogenDependentCellCycleModel::SetMorphogenInfluence(double morphogenInfluence)
{
	mMorphogenInfluence = morphogenInfluence;
}

double MorphogenDependentCellCycleModel::GetMorphogenInfluence()
{
    return mMorphogenInfluence;
}


void MorphogenDependentCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<TargetMass>" << mTargetMass << "</TargetMass>\n";
    *rParamsFile << "\t\t\t<MorphogenInfluence>" << mMorphogenInfluence << "</MorphogenInfluence>\n";

    // Call method on direct parent class
    AbstractCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(MorphogenDependentCellCycleModel)
