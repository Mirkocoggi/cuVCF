#!/bin/bash

#/usr/bin/time -f "Max Mem: %M KB" ./bin/VCFparser -v data/tiny2.vcf -t 1
#
# Define the number of times to run the C++ program
NUM_RUNS=3

NUM_TH=24
make GPU
# Loop to run the C++ program multiple times
echo "IRBT3M:"
echo "IRBT3M:" >> GenoGraMachine.txt

for((j = 1; j <= NUM_TH; j=j*2)); do
    echo "Th $j:"
    echo "Th $j:" >> GenoGraMachine.txt
    for ((i = 1; i <= NUM_RUNS; i++)); do
        echo "Run $i"   
        /usr/bin/time -f "Time: %E | Max Mem: %M KB" ./bin/VCFparser -v data/IRBT3M.vcf -t $j >> GenoGraMachine.txt 2>&1
        #echo "-----" >>  GenoGraMachine.txt
    done

    echo "===================\n" >>  GenoGraMachine.txt
done

echo "IRBT:"
echo "IRBT:" >> GenoGraMachine.txt

for((j = 1; j <= NUM_TH; j=j*2)); do
    echo "Th $j:"
    echo "Th $j:" >> GenoGraMachine.txt
    for ((i = 1; i <= NUM_RUNS; i++)); do
        echo "Run $i"   
        /usr/bin/time -f "Time: %E | Max Mem: %M KB" ./bin/VCFparser -v data/IRBT.vcf -t $j >> GenoGraMachine.txt 2>&1
        #echo "-----" >>  GenoGraMachine.txt
    done

    echo "===================\n" >>  GenoGraMachine.txt
done

echo "danio_rerio:"
echo "danio_rerio:" >> GenoGraMachine.txt

for((j = 1; j <= NUM_TH; j=j*2)); do
    echo "Th $j:"
    echo "Th $j:" >> GenoGraMachine.txt
    for ((i = 1; i <= NUM_RUNS; i++)); do
        echo "Run $i"   
        /usr/bin/time -f "Time: %E | Max Mem: %M KB" ./bin/VCFparser -v data/danio_rerio.vcf -t $j >> GenoGraMachine.txt 2>&1
        #echo "-----" >>  GenoGraMachine.txt
    done

    echo "===================\n" >>  GenoGraMachine.txt
done

echo "felis_catus:"
echo "felis_catus:" >> GenoGraMachine.txt

for((j = 1; j <= NUM_TH; j=j*2)); do
    echo "Th $j:"
    echo "Th $j:" >> GenoGraMachine.txt
    for ((i = 1; i <= NUM_RUNS; i++)); do
        echo "Run $i"   
        /usr/bin/time -f "Time: %E | Max Mem: %M KB" ./bin/VCFparser -v data/felis_catus.vcf -t $j >> GenoGraMachine.txt 2>&1
        #echo "-----" >>  GenoGraMachine.txt
    done

    echo "===================\n" >>  GenoGraMachine.txt
done

echo "bos_taurus:"
echo "bos_taurus:" >> GenoGraMachine.txt

for((j = 1; j <= NUM_TH; j=j*2)); do
    echo "Th $j:"
    echo "Th $j:" >> GenoGraMachine.txt
    for ((i = 1; i <= NUM_RUNS; i++)); do
        echo "Run $i"   
        /usr/bin/time -f "Time: %E | Max Mem: %M KB" ./bin/VCFparser -v data/bos/bos_taurus.vcf -t $j >> GenoGraMachine.txt 2>&1
        #echo "-----" >>  GenoGraMachine.txt
    done

    echo "===================\n" >>  GenoGraMachine.txt
done


echo "bos75M:"
echo "bos75M:" >> GenoGraMachine.txt

for((j = 1; j <= NUM_TH; j=j*2)); do
    echo "Th $j:"
    echo "Th $j:" >> GenoGraMachine.txt
    for ((i = 1; i <= NUM_RUNS; i++)); do
        echo "Run $i"   
        /usr/bin/time -f "Time: %E | Max Mem: %M KB" ./bin/VCFparser -v data/bos75M.vcf -t $j >> GenoGraMachine.txt 2>&1
        #echo "-----" >>  GenoGraMachine.txt
    done

    echo "===================\n" >>  GenoGraMachine.txt
done

#homo_sapiens-chr20.vcf