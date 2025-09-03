#!/bin/bash

#/usr/bin/time -f "Max Mem: %M KB" ./bin/VCFparser -v data/tiny2.vcf -t 1
#
# Define the number of times to run the C++ program
NUM_RUNS=5

NUM_TH=24
make GPU
# Loop to run the C++ program multiple times
echo "IRBT3M:"
echo "IRBT3M:" >> output.txt

for((j = 1; j <= NUM_TH; j=j*2)); do
    echo "Th $j:"
    echo "Th $j:" >> output.txt
    for ((i = 1; i <= NUM_RUNS; i++)); do
        echo "Run $i"   
        /usr/bin/time -f "Time: %E | Max Mem: %M KB" ./bin/VCFparser -v data/IRBT3M.vcf -t $j >> output.txt 2>&1
        
    done

    echo "===================\n" >>  output.txt
done

echo "IRBT:"
echo "IRBT:" >> output.txt

for((j = 1; j <= NUM_TH; j=j*2)); do
    echo "Th $j:"
    echo "Th $j:" >> output.txt
    for ((i = 1; i <= NUM_RUNS; i++)); do
        echo "Run $i"   
        /usr/bin/time -f "Time: %E | Max Mem: %M KB" ./bin/VCFparser -v data/IRBT.vcf -t $j >> output.txt 2>&1
        
    done

    echo "===================\n" >>  output.txt
done

echo "danio_rerio:"
echo "danio_rerio:" >> output.txt

for((j = 1; j <= NUM_TH; j=j*2)); do
    echo "Th $j:"
    echo "Th $j:" >> output.txt
    for ((i = 1; i <= NUM_RUNS; i++)); do
        echo "Run $i"   
        /usr/bin/time -f "Time: %E | Max Mem: %M KB" ./bin/VCFparser -v data/danio_rerio.vcf -t $j >> output.txt 2>&1
        
    done

    echo "===================\n" >>  output.txt
done

echo "felis_catus:"
echo "felis_catus:" >> output.txt

for((j = 1; j <= NUM_TH; j=j*2)); do
    echo "Th $j:"
    echo "Th $j:" >> output.txt
    for ((i = 1; i <= NUM_RUNS; i++)); do
        echo "Run $i"   
        /usr/bin/time -f "Time: %E | Max Mem: %M KB" ./bin/VCFparser -v data/felis_catus.vcf -t $j >> output.txt 2>&1
        
    done

    echo "===================\n" >>  output.txt
done

echo "bos_taurus:"
echo "bos_taurus:" >> output.txt

for((j = 1; j <= NUM_TH; j=j*2)); do
    echo "Th $j:"
    echo "Th $j:" >> output.txt
    for ((i = 1; i <= NUM_RUNS; i++)); do
        echo "Run $i"   
        /usr/bin/time -f "Time: %E | Max Mem: %M KB" ./bin/VCFparser -v data/bos/bos_taurus.vcf -t $j >> output.txt 2>&1
        
    done

    echo "===================\n" >>  output.txt
done


echo "bos75M:"
echo "bos75M:" >> output.txt

for((j = 1; j <= NUM_TH; j=j*2)); do
    echo "Th $j:"
    echo "Th $j:" >> output.txt
    for ((i = 1; i <= NUM_RUNS; i++)); do
        echo "Run $i"   
        /usr/bin/time -f "Time: %E | Max Mem: %M KB" ./bin/VCFparser -v data/bos75M.vcf -t $j >> output.txt 2>&1
        
    done

    echo "===================\n" >>  output.txt
done
