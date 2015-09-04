#!/bin/bash
#
# Script to compile an executable and then run a sweep by calling the run script
#

cd ../../../

scons co=1 b=GccOpt ts=projects/CellBasedComparison/test/TestCellSorting.hpp

cd projects/CellBasedComparison/test


echo "Build Complete"

/bin/bash run_sorting.sh