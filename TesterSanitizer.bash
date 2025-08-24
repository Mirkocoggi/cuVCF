#!/bin/bash

#make DEBUG;
echo "run sanitizer"
#compute-sanitizer ./bin/VCFparser -v data/bos_head.vcf -t 8 > h.log
#compute-sanitizer ./bin/VCFparser -v data/IRBT.vcf -t 16 > h.log
cuda-gdb --args ./bin/VCFparser -v data/chrx_AAAAAA.vcf -t 1
#echo $? >> h.log

#If needed use cuda-gdb as follow to find segfault and backtrace the error:
#cuda-gdb --args ./bin/VCFparser -v data/IRBT.vcf -t 4
#then: "run" then: "bt" to backtrace