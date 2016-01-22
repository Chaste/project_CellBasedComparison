#!/bin/bash
#
# Script to illustrate running batch jobs and passing in arguments.
#
# 
# This script assumes that the following has been run successfully:
# scons co=1 b=GccOpt ts=projects/Ozzy/test/CellBasedComparison/TestCylindricalCrypt.hpp
#

num_sims=1;

#CI_LEVEL[0]="1.0"
#CI_LEVEL[0]="0.9"
CI_LEVEL[0]="0.8"
#CI_LEVEL[2]="0.7"
#CI_LEVEL[4]="0.6"
#CI_LEVEL[5]="0.5"

for (( i=0 ; i<${num_sims} ; i++))
do
    echo "Run " $i;
    for (( j=0 ; j<${#CI_LEVEL[*]} ; j++))
    do
    	echo "  CI level " ${CI_LEVEL[$j]};
    	# NB "nice -20" gives the jobs low priority (good if they are going to dominate the server and no slower if nothing else is going on)
    	# ">" directs std::cout to the file.
    	# "2>&1" directs std::cerr to the same place.
    	# "&" on the end lets the script carry on and not wait until this has finished.
    	nice -20 ../build/optimised/TestCylindricalCryptRunner -sim_index $i -CI ${CI_LEVEL[$j]} > output/CryptRun_${i}_${CI_LEVEL[$j]}_Output.txt 2>&1 &
    done
done

echo "Jobs submitted"