#!/bin/bash

#SBATCH -n 8                               # 1 core
#SBATCH -t 0-05:00                         # Runtime of 5 minutes, in D-HH:MM format
#SBATCH --mem=30G
#SBATCH -p park                           # Run in short partition
#SBATCH -o hostname_%j.out                 # File to which STDOUT + STDERR will be written, including job ID in filename
#SBATCH --mail-type=END                    # ALL email notification type
#SBATCH --mail-user=chong.simonchu@gmail.com  # Email to which notifications will be sent
#SBATCH --account=park_contrib
#
##NA12878_illumina.list ./NA12878/  8 ./NA12878/CEUTrio.HiSeq.WGS.b37_decoy.NA12878.clean.dedup.recal.20120117-ram/CEUTrio.HiSeq.WGS.b37_decoy.NA12878.clean.dedup.recal.20120117.tea.L1   /n/data1/hms/dbmi/park/simon_chu/projects/tools/MELTv2.1.4/me_refs/LINE1_MELT/LINE1.fa /n/data1/hms/dbmi/park/SOFTWARE/LongRanger/refdata-b37-2.1.0/fasta/genome.fa NA12878_illumina.cns.list
#
BAM_LIST=NA12878_illumina.list
TMP_CNS=./NA12878/
CORES=8
INPUT=./NA12878/CEUTrio.HiSeq.WGS.b37_decoy.NA12878.clean.dedup.recal.20120117-ram/CEUTrio.HiSeq.WGS.b37_decoy.NA12878.clean.dedup.recal.20120117.tea.L1
L1_CNS=/n/data1/hms/dbmi/park/simon_chu/projects/tools/MELTv2.1.4/me_refs/LINE1_MELT/LINE1.fa
REF=/n/data1/hms/dbmi/park/SOFTWARE/LongRanger/refdata-b37-2.1.0/fasta/genome.fa
OUTPUT=NA12878_illumina.cns.list
#
python ./../filter_src/TEA_filter.py -b ${BAM_LIST} -p ${TMP_CNS} -n ${CORES} -i ${INPUT}  -c ${L1_CNS} --r ${REF}  -o ${OUTPUT}
