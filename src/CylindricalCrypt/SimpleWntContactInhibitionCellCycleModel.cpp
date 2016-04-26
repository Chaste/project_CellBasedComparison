/*

Copyright (c) 2005-2015, University of Oxford.
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

#include "SimpleWntContactInhibitionCellCycleModel.hpp"
#include "CellLabel.hpp"
#include "RandomNumberGenerator.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "WntConcentration.hpp"

SimpleWntContactInhibitionCellCycleModel::SimpleWntContactInhibitionCellCycleModel()
    : ContactInhibitionCellCycleModel(),
      mWntThreshold(0.5)
{
}

SimpleWntContactInhibitionCellCycleModel::SimpleWntContactInhibitionCellCycleModel(const SimpleWntContactInhibitionCellCycleModel& rModel)
   : ContactInhibitionCellCycleModel(rModel),
     mWntThreshold(rModel.mWntThreshold)
{
    /*
     * Set each member variable of the new cell-cycle model that inherits
     * its value from the parent.
     *
     * Note 1: some of the new cell-cycle model's member variables will already
     * have been correctly initialized in its constructor or parent classes.
     *
     * Note 2: one or more of the new cell-cycle model's member variables
     * may be set/overwritten as soon as InitialiseDaughterCell() is called on
     * the new cell-cycle model.
     *
     * Note 3: Only set the variables defined in this class. Variables defined
     * in parent classes will be defined there.
     *
     */
}

void SimpleWntContactInhibitionCellCycleModel::UpdateCellCyclePhase()
{
    double wnt_level= GetWntLevel();

    // Set the cell type to TransitCellProliferativeType if the Wnt stimulus exceeds wnt_division_threshold if not set it to Differentiated Type
    if (wnt_level >= GetWntThreshold())
    {
        boost::shared_ptr<AbstractCellProperty> p_transit_type =
            mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<TransitCellProliferativeType>();
        mpCell->SetCellProliferativeType(p_transit_type);
    }
    else
    {
        // The cell is set to have DifferentiatedCellProliferativeType and so in G0 phase
        boost::shared_ptr<AbstractCellProperty> p_diff_type =
            mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<DifferentiatedCellProliferativeType>();
        mpCell->SetCellProliferativeType(p_diff_type);
    }

    ContactInhibitionCellCycleModel::UpdateCellCyclePhase();
}

void SimpleWntContactInhibitionCellCycleModel::SetG1Duration()
{
    assert(mpCell != NULL);

    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();

    if (mpCell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>())
    {
        mG1Duration = p_gen->NormalRandomDeviate(GetTransitCellG1Duration(), 1.0);
    }
    else if (mpCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())
    {
        mG1Duration = DBL_MAX;
    }
    else
    {
        NEVER_REACHED;
    }

    // Check that the normal random deviate has not returned a small or negative G1 duration
    if (mG1Duration < mMinimumGapDuration)
    {
        mG1Duration = mMinimumGapDuration;
    }
}

double SimpleWntContactInhibitionCellCycleModel::GetWntLevel()
{
    assert(mpCell != NULL);
    double level = 0;

    switch (mDimension)
    {
        case 1:
        {
            const unsigned DIM = 1;
            level = WntConcentration<DIM>::Instance()->GetWntLevel(mpCell);
            break;
        }
        case 2:
        {
            const unsigned DIM = 2;
            level = WntConcentration<DIM>::Instance()->GetWntLevel(mpCell);
            break;
        }
        case 3:
        {
            const unsigned DIM = 3;
            level = WntConcentration<DIM>::Instance()->GetWntLevel(mpCell);
            break;
        }
        default:
            NEVER_REACHED;
    }
    return level;
}

AbstractCellCycleModel* SimpleWntContactInhibitionCellCycleModel::CreateCellCycleModel()
{
    // Create a new cell-cycle model
    return new SimpleWntContactInhibitionCellCycleModel(*this);
}

void SimpleWntContactInhibitionCellCycleModel::SetWntThreshold(double wntThreshold)
{
    mWntThreshold = wntThreshold;
}

double SimpleWntContactInhibitionCellCycleModel::GetWntThreshold()
{
    return mWntThreshold;
}

void SimpleWntContactInhibitionCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<WntThreshold>" << mWntThreshold << "</WntThreshold>\n";

    // Call method on direct parent class
    ContactInhibitionCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(SimpleWntContactInhibitionCellCycleModel)
