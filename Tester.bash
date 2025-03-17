#!/bin/bash

#/usr/bin/time -f "Max Mem: %M KB" ./bin/VCFparser -v data/tiny2.vcf -t 1
#
# Define the number of times to run the C++ program
NUM_RUNS=3

NUM_TH=24

# Loop to run the C++ program multiple times
for((j = 1; j <= NUM_TH; j=j*2)); do
    echo "Th $j:"
    echo "Th $j:" >> GenoGraMachine.txt
    for ((i = 1; i <= NUM_RUNS; i++)); do
        echo "Run $i"   
        /usr/bin/time -f "Time: %E | Max Mem: %M KB" ./bin/VCFparser -v data/IRBT3M.vcf -t $j >> GenoGraMachineSamples.txt 2>&1
        #echo "-----" >>  GenoGraMachine.txt
    done

    echo "===================\n" >>  GenoGraMachine.txt
done

#homo_sapiens-chr20.vcf