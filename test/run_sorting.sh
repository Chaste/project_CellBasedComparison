#!/bin/bash
#
# Script to illustrate running batch jobs and passing in arguments.
#
# 
# This script assumes that the following has been run successfully:
# scons co=1 b=GccOpt ts=projects/Ozzy/test/CellBasedComparison/TestCylindricalCrypt.hpp
#

start_sim=9;
end_sim=10;

CCD[0]="10.0"
CCD[1]="50.0"
CCD[2]="100.0"

for (( i=${start_sim} ; i<${end_sim} ; i++))
do
    echo "Run " $i;
    for (( j=0 ; j<${#CCD[*]} ; j++))
    do
    	echo " Mean  CCD " ${CCD[$j]};
    	# NB "nice -20" gives the jobs low priority (good if they are going to dominate the server and no slower if nothing else is going on)
    	# ">" directs std::cout to the file.
    	# "2>&1" directs std::cerr to the same place.
    	# "&" on the end lets the script carry on and not wait until this has finished.
    	nice -20 ../build/optimised/TestCellSortingRunner -sim_index $i -CCD ${CCD[$j]} > output/SortingRun_${i}_${CCD[$j]}_Output.txt 2>&1 &
    done
done

echo "Jobs submitted"