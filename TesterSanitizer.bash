#!/bin/bash

#cmake -S . -B build-dbg -DCMAKE_BUILD_TYPE=Debug -DCUVCF_CUDA=ON && cmake --build build-dbg --target VCFparser_gpu -j
echo "run sanitizer"
#compute-sanitizer ./build-dbg/VCFparser_gpu -v data/bos_head.vcf -t 8 > h.log
#compute-sanitizer ./build-dbg/VCFparser_gpu -v data/IRBT.vcf -t 16 > h.log
cuda-gdb --args ./build-dbg/VCFparser_gpu -v data/chrx_AAAAAA.vcf -t 1
#echo $? >> h.log

#If needed use cuda-gdb as follow to find segfault and backtrace the error:
#cuda-gdb --args ./build-dbg/VCFparser_gpu -v data/IRBT.vcf -t 4
#then: "run" then: "bt" to backtrace