#!/bin/bash

./TestFiltriNoSamp_FELIS.bash > out.log 2>&1
./TestFiltriNoSamp_DANIO.bash >> out.log 2>&1
./TestFiltriNoSamp_BOS.bash >> out.log 2>&1
./TesterFiltriSamp.bash >> out.log 2>&1