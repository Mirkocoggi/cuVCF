#!/bin/bash

# Define the number of times to run the C++ program
NUM_RUNS=30
NUM_TH=64

# Loop to run the C++ program multiple times
for((j = 1; j <= NUM_TH; j=j*2)); do

    echo "Th $j:" >> ResPopulateCR20.txt
    for ((i = 1; i <= NUM_RUNS; i++)); do
        #echo "Run $i:" >> results.txt
        ./bin/VCFparser -v data/homo_sapiens-chr20.vcf -t $j >>  ResPopulateCR20.txt
    done

    echo "===================" >>  ResPopulateCR20.txt
done

#homo_sapiens-chr20.vcf