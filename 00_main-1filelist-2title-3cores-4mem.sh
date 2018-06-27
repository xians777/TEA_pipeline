#!/bin/bash

FILELIST=$1
GROUP_NAME=$2
CORES=$3
MEM=$4

cat ${FILELIST} | xargs -I {} ./01_sbatch_run_tea-1file-2title-3cores-4mem.sh {} ${GROUP_NAME} ${CORES} ${MEM}
