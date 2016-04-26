#!/bin/bash
#
# Script to compile an executable and then run a sweep by calling the run script
#

cd ../../../../

scons co=1 b=GccOpt -j6 ts=projects/CellBasedComparison/test/Sweeps/TestMorphogenMonolayerSweeps.hpp

cd projects/CellBasedComparison/test/Sweeps


echo "Build Complete"

/bin/bash run_morphogen.sh