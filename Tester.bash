#!/bin/bash

# Define the number of times to run the C++ program
NUM_RUNS=30
NUM_TH=12

# Loop to run the C++ program multiple times
for((j = 1; j <= NUM_TH; j++)); do

    echo "Th $j:" >> resultsArr.txt

    for ((i = 1; i <= NUM_RUNS; i++)); do
        #echo "Run $i:" >> results.txt
        ./bin/VCFparser -v data/test_1600.vcf -t $j >> resultsArr.txt
    done

    echo "===================" >> resultsArr.txt
done