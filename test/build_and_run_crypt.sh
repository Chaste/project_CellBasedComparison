#!/bin/bash
#
# Script to compile an executable and then run a sweep by calling the run script
#

cd ../../../

scons co=1 b=GccOpt -j6 ts=projects/CellBasedComparison/test/TestCylindricalCrypt.hpp

cd projects/CellBasedComparison/test


echo "Build Complete"

/bin/bash run_crypt.sh