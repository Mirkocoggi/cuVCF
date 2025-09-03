#!/bin/bash

./TestFilterNoSamp_FELIS.bash > out.log 2>&1
./TestFilterNoSamp_DANIO.bash >> out.log 2>&1
./TestFilterNoSamp_BOS.bash >> out.log 2>&1
./TesterFilterSamp.bash >> out.log 2>&1