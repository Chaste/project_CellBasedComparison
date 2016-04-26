#!/bin/bash
#
# Script to illustrate running batch jobs and passing in arguments.
#
# 
# This script assumes that the following has been run successfully:
# scons co=1 b=GccOpt ts=projects/Ozzy/test/CellBasedComparison/TestDeltaNotch.hpp
#

num_sims=10;

PDIV[0]="0.1"
PDIV[1]="0.05";
PDIV[2]="0.01";

for (( i=0 ; i<${num_sims} ; i++))
do
    echo "Run " $i;
    for (( j=0 ; j<${#PDIV[*]} ; j++))
    do
    	echo " prob division " ${PDIV[$j]};
    	# NB "nice -20" gives the jobs low priority (good if they are going to dominate the server and no slower if nothing else is going on)
    	# ">" directs std::cout to the file.
    	# "2>&1" directs std::cerr to the same place.
    	# "&" on the end lets the script carry on and not wait until this has finished.
    	nice -20 ../../build/optimised/Sweeps/TestDeltaNotchSweepsRunner -sim_index $i -prob_division ${PDIV[$j]} > output/DeltaRun_${i}_${PDIV[$j]}_Output.txt 2>&1 &
    done
done

echo "Jobs submitted"