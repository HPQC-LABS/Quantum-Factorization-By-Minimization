#!/bin/sh

cd /home/tanburn/Quantum-Factorization-By-Minimization/batch_scripts_hamming


if [ "$1" == "" ]; then
    exit 1
fi

EXP_NAME="hexp$1"

screen -S $EXP_NAME -d -m

BATCH_FILE="batch_experiments_$1x$1.sh"

screen -S $EXP_NAME -X stuff ". ./$BATCH_FILE\r"