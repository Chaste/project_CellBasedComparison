#!/bin/bash
#
# Script to illustrate running batch jobs and passing in arguments.
#
# 
# This script assumes that the following has been run successfully:
# scons co=1 b=GccOpt ts=projects/Ozzy/test/CellBasedComparison/TestCylindricalCrypt.hpp
#

start_sim=0;
end_sim=1;

noise[0]="100"
noise[1]="10"
noise[2]="1"
noise[3]="10"
noise[4]="100"


for (( i=${start_sim} ; i<${end_sim} ; i++))
do
    echo "Run " $i;
    for (( j=0 ; j<${#noise[*]} ; j++))
    do
    	echo " Mean  noise " ${noise[$j]};
    	# NB "nice -20" gives the jobs low priority (good if they are going to dominate the server and no slower if nothing else is going on)
    	# ">" directs std::cout to the file.
    	# "2>&1" directs std::cerr to the same place.
    	# "&" on the end lets the script carry on and not wait until this has finished.
    	nice -20 ../build/optimised/TestCellSortingRunner -sim_index $i -noise ${noise[$j]} > output/SortingRun_${i}_${noise[$j]}_Output.txt 2>&1 &
    done
done

echo "Jobs submitted"