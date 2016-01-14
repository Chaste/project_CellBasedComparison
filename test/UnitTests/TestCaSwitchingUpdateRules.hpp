/*

Copyright (c) 2005-2016, University of Oxford.
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
#ifndef TESTCASWITCHINGUPDATERULES_HPP_
#define TESTCASWITCHINGUPDATERULES_HPP_

#include <cxxtest/TestSuite.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>


#include "AbstractCaSwitchingUpdateRule.hpp"
#include "DifferentialAdhesionCaSwitchingUpdateRule.hpp"
#include "DifferentialAdhesionCaSwitchingUpdateRule.hpp"
#include "CellsGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "CaBasedCellPopulation.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "WildTypeCellMutationState.hpp"
#include "CellLabel.hpp"
#include "PottsMeshGenerator.hpp"
#include "NodesOnlyMesh.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "SmartPointers.hpp"
#include "FileComparison.hpp"

#include "PetscSetupAndFinalize.hpp"

class TestCaSwitchingUpdateRules : public AbstractCellBasedTestSuite
{
public:
    /*
     * Now test the switching rules.
     */
    void TestDifferentialAdhesionCaSwitchingUpdateRuleIn2d() throw (Exception)
    {
        // Set the timestep and size of domain to let us calculate the probabilities of movement
        double delta_t = 0.1;
        double delta_x = 1;

        // Create an update law system
        DifferentialAdhesionCaSwitchingUpdateRule<2> random_switching_update_rule;

        // Test get/set methods
        TS_ASSERT_DELTA(random_switching_update_rule.GetLabelledCellLabelledCellAdhesionEnergyParameter(), 0.2, 1e-12);
        TS_ASSERT_DELTA(random_switching_update_rule.GetLabelledCellCellAdhesionEnergyParameter(), 0.1, 1e-12);
        TS_ASSERT_DELTA(random_switching_update_rule.GetCellCellAdhesionEnergyParameter(), 0.2, 1e-12);
        TS_ASSERT_DELTA(random_switching_update_rule.GetCellBoundaryAdhesionEnergyParameter(), 0.2, 1e-12);
        TS_ASSERT_DELTA(random_switching_update_rule.GetLabelledCellBoundaryAdhesionEnergyParameter(), 0.2, 1e-12);
        TS_ASSERT_DELTA(random_switching_update_rule.GetTemperature(), 0.1, 1e-12);


        random_switching_update_rule.SetLabelledCellLabelledCellAdhesionEnergyParameter(1.0);
        random_switching_update_rule.SetLabelledCellCellAdhesionEnergyParameter(2.0);
        random_switching_update_rule.SetCellCellAdhesionEnergyParameter(3.0);
        random_switching_update_rule.SetCellBoundaryAdhesionEnergyParameter(4.0);
        random_switching_update_rule.SetLabelledCellBoundaryAdhesionEnergyParameter(5.0);
        random_switching_update_rule.SetTemperature(10);


        TS_ASSERT_DELTA(random_switching_update_rule.GetLabelledCellLabelledCellAdhesionEnergyParameter(), 1.0, 1e-12);
        TS_ASSERT_DELTA(random_switching_update_rule.GetLabelledCellCellAdhesionEnergyParameter(), 2.0, 1e-12);
        TS_ASSERT_DELTA(random_switching_update_rule.GetCellCellAdhesionEnergyParameter(), 3.0, 1e-12);
        TS_ASSERT_DELTA(random_switching_update_rule.GetCellBoundaryAdhesionEnergyParameter(), 4.0, 1e-12);
        TS_ASSERT_DELTA(random_switching_update_rule.GetLabelledCellBoundaryAdhesionEnergyParameter(), 5.0, 1e-12);
        TS_ASSERT_DELTA(random_switching_update_rule.GetTemperature(), 10.0, 1e-12);

        // Test EvaluateProbability()

        // Create a simple 2D PottsMesh
        PottsMeshGenerator<2> generator(4, 0, 0, 4, 0, 0);
        PottsMesh<2>* p_mesh = generator.GetMesh();

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<FixedDurationGenerationBasedCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, 4u, p_diff_type);

        // Specify where cells lie here we have 2 cells in the centre on each of the on the bottom two rows
        // the left two are labeled
        std::vector<unsigned> location_indices;

        location_indices.push_back(1);
        location_indices.push_back(2);

        location_indices.push_back(5);
        location_indices.push_back(6);

        // Label the cells 0 and 2
        boost::shared_ptr<AbstractCellProperty> p_label(CellPropertyRegistry::Instance()->Get<CellLabel>());
        cells[0]->AddCellProperty(p_label);
        cells[2]->AddCellProperty(p_label);

        // Create cell population
        CaBasedCellPopulation<2u> cell_population(*p_mesh, cells, location_indices);

        // Labeled cells with empty nodes
        TS_ASSERT_DELTA(random_switching_update_rule.EvaluateHamiltonian(0,1,cell_population),11.0,1e-6);
        TS_ASSERT_DELTA(random_switching_update_rule.EvaluateSwitchingProbability(0,1,cell_population,delta_t,delta_x),delta_t*exp(-11.0/10),1e-6);
        TS_ASSERT_DELTA(random_switching_update_rule.EvaluateHamiltonian(4,5,cell_population),11,1e-6);
        TS_ASSERT_DELTA(random_switching_update_rule.EvaluateHamiltonian(5,8,cell_population),11,1e-6);
        TS_ASSERT_DELTA(random_switching_update_rule.EvaluateHamiltonian(5,9,cell_population),16,1e-6);

        // Non Labeled cells with empty nodes
        TS_ASSERT_DELTA(random_switching_update_rule.EvaluateHamiltonian(2,3,cell_population),8,1e-6);
        TS_ASSERT_DELTA(random_switching_update_rule.EvaluateHamiltonian(6,7,cell_population),8,1e-6);
        TS_ASSERT_DELTA(random_switching_update_rule.EvaluateHamiltonian(6,10,cell_population),12,1e-6);

        // Same cell types
        TS_ASSERT_DELTA(random_switching_update_rule.EvaluateHamiltonian(1,5,cell_population),0,1e-6);
        TS_ASSERT_DELTA(random_switching_update_rule.EvaluateHamiltonian(2,6,cell_population),0,1e-6);
        TS_ASSERT_DELTA(random_switching_update_rule.EvaluateSwitchingProbability(2,6,cell_population,delta_t,delta_x),delta_t,1e-6);


        // Labeled Cells with Non labeled ones
        TS_ASSERT_DELTA(random_switching_update_rule.EvaluateHamiltonian(5,6,cell_population),0,1e-6);// 0 as to 3+1 = 2 + 2
        TS_ASSERT_DELTA(random_switching_update_rule.EvaluateHamiltonian(1,2,cell_population),0,1e-6);// 0 as to 3+1 = 2 + 2
        TS_ASSERT_DELTA(random_switching_update_rule.EvaluateHamiltonian(1,6,cell_population),1,1e-6);
        TS_ASSERT_DELTA(random_switching_update_rule.EvaluateHamiltonian(2,5,cell_population),-1,1e-6);
        TS_ASSERT_DELTA(random_switching_update_rule.EvaluateSwitchingProbability(2,5,cell_population,delta_t,delta_x),delta_t,1e-6);



    }

    void noTestArchiveDifferentialAdhesionCaSwitchingUpdateRule() throw(Exception)
    {
        OutputFileHandler handler("archive", false);
        std::string archive_filename = handler.GetOutputDirectoryFullPath() + "DifferentialAdhesionCaSwitchingUpdateRule.arch";

        {
            DifferentialAdhesionCaSwitchingUpdateRule<2> update_rule;

            std::ofstream ofs(archive_filename.c_str());
            boost::archive::text_oarchive output_arch(ofs);

            // Set member variables
            //update_rule.SetSwitchingParameter(1.0);

            // Serialize via pointer to most abstract class possible
            AbstractCaSwitchingUpdateRule<2>* const p_update_rule = &update_rule;
            output_arch << p_update_rule;
        }

        {
            AbstractCaSwitchingUpdateRule<2>* p_update_rule;

            // Create an input archive
            std::ifstream ifs(archive_filename.c_str(), std::ios::binary);
            boost::archive::text_iarchive input_arch(ifs);

            // Restore from the archive
            input_arch >> p_update_rule;

            // Test the member data
            //TS_ASSERT_DELTA((static_cast<DifferentialAdhesionCaSwitchingUpdateRule<2>*>(p_update_rule))->GetSwitchingParameter(), 1.0, 1e-6);

            // Tidy up
            delete p_update_rule;
        }
    }

    void noTestSwitchingUpdateRuleOutputUpdateRuleInfo()
    {
        std::string output_directory = "TestCaSwitchingUpdateRulesOutputParameters";
        OutputFileHandler output_file_handler(output_directory, false);

        // Test with DifferentialAdhesionCaSwitchingUpdateRule
        DifferentialAdhesionCaSwitchingUpdateRule<2> random_switching_update_rule;
        //random_switching_update_rule.SetSwitchingParameter(1.0);

        TS_ASSERT_EQUALS(random_switching_update_rule.GetIdentifier(), "DifferentialAdhesionCaSwitchingUpdateRule-2");

        out_stream random_switching_update_rule_parameter_file = output_file_handler.OpenOutputFile("random_switching_update_rule_results.parameters");
        random_switching_update_rule.OutputUpdateRuleInfo(random_switching_update_rule_parameter_file);
        random_switching_update_rule_parameter_file->close();

        // Compare the generated file in test output with a reference copy in the source code.
        FileFinder generated = output_file_handler.FindFile("random_switching_update_rule_results.parameters");
        FileFinder reference("cell_based/test/data/TestCaUpdateRules/random_switching_update_rule_results.parameters",
                RelativeTo::ChasteSourceRoot);
        FileComparison comparer(generated, reference);
        TS_ASSERT(comparer.CompareFiles());
    }
};

#endif /*TESTCASWITCHINGUPDATERULES_HPP_*/
