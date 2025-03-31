#!/bin/bash

make DEBUG;
echo "run sanitizer"
#compute-sanitizer ./bin/VCFparser -v data/bos_head.vcf -t 8 > h.log
compute-sanitizer ./bin/VCFparser -v data/IRBT3M.vcf -t 16 > h.log
echo $? >> h.log

#If needed use cuda-gdb as follow to find segfault and backtrace the error:
#cuda-gdb --args ./bin/VCFparser -v data/bos75M.vcf -t 8
#then: "run" then: "bt" to backtrace