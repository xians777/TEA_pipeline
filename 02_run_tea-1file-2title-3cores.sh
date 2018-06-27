#!/bin/bash -x

module load gcc/6.2.0
module load boost/1.62.0
module load pigz/2.3.4
module load R/3.4.1

module load samtools/1.3.1
module load bedtools/2.26.0

# for sambamba and bwa
for REQ in /n/data1/bch/genetics/lee/kyu/opt/bin/ /n/app/bwa/0.6.2/bin/
do
	if [[ ${PATH} != *"${REQ}"* ]];then
    		export PATH=${REQ}:$PATH
	fi
done

# gcc and mybamtools
for REQ in /n/app/gcc/6.2.0/lib64 /n/app/gcc/6.2.0/lib /home/ly55/opt/meerkat/lib
do
	if [[ ${LD_LIBRARY_PATH} != *"${REQ}"* ]];then
		export LD_LIBRARY_PATH=${REQ}:$LD_LIBRARY_PATH
	fi
done

export PERL5LIB=/home/ly55/perl/share/perl/5.10.0:/home/ly55/perl/lib/perl5/x86_64-linux-gnu-thread-multi:/home/ly55/perl/lib/perl/5.10.0:/home/el114/perl5/BioPerl-1.6.924

# this will be read by TEA
export tea_base=/n/data1/bch/genetics/lee/kyu/opt/tea

# check if bamsize is not zero
if [ -e ${1} ]; then
	BAM_SIZE=$(stat -Lc%s "${1}")
    if (( ${BAM_SIZE} < 1 )); then
    	exit 0
    fi
fi

# default Group Directory name is "tea"
if [ -z $2 ]; then
	DIR_GROUP="tea"
else
	DIR_GROUP=$2
fi

# default number of core is 4
if [ -z $3 ]; then
	CORES=4
else
	CORES=$3
fi

BAM_FILE=$1
DIR_NAME=$( dirname ${BAM_FILE} )
BASE_NAME=$( basename ${BAM_FILE} )
PREFIX=${BASE_NAME%.bam}
INPUT_PREFIX=${1%.bam}

LOGS_DIR=${DIR_GROUP}/logs
LOGS_FILE=${LOGS_DIR}/${PREFIX}

mkdir -p ${LOGS_DIR}
mkdir -p ${DIR_GROUP}/${PREFIX}/${PREFIX}


NAME_REF="hs37d5"
DIR_REF=/n/data1/bch/genetics/lee/data/references/bwa-0.6.2-r126
PATH_REF=${DIR_REF}/${NAME_REF}
PATH_REF_FASTA=${DIR_REF}/${NAME_REF}/${NAME_REF}.fasta
#PATH_REF_FASTA=/n/data1/bch/genetics/lee/Rebeca/Tea_Melt_benchmark_6_6_18/reference/genome.fa
PATH_REF_FAI=${PATH_REF_FASTA}.fai

PATH_REFERENCE_FA=/n/data1/bch/genetics/lee/kyu/opt/vince/tea/lib/assembly/human_youngTE_revisedPolyA.fa
PATH_RMASKER_TXT=/n/data1/bch/genetics/lee/kyu/opt/vince/tea/lib/rmasker/rmasker.hg19.youngte.RData.txt

PATH_MEERKAT=/n/data1/bch/genetics/lee/kyu/opt/meerkat

time perl ${PATH_MEERKAT}/scripts/pre_process.pl  \
    -D ${DIR_GROUP} -b ${BAM_FILE} -I ${PATH_REF_FASTA} -A ${PATH_REF_FAI} -t ${CORES} \
    -u 1 -s 20 -k 1500 -q 15 \
    >${LOGS_FILE}.job1log 2>${LOGS_FILE}.job1err

time perl ${PATH_MEERKAT}/scripts/meerkat.pl -P dc  \
    -D ${DIR_GROUP} -b ${BAM_FILE} -F ${PATH_REF_FASTA} -t ${CORES} \
    -s 20 -d 5 -p 3 -o 1 -m 0  \
    >${LOGS_FILE}.job2log 2>${LOGS_FILE}.job2err

TEA_EXE=/n/data1/bch/genetics/lee/kyu/opt/meerkat/bin/tea

time ${TEA_EXE} \
    tea \
    -D ${DIR_GROUP} -b ${BAM_FILE} -r ${PATH_REFERENCE_FA} -M ${PATH_RMASKER_TXT} \
    --mem --force --include_head_clip --oi --contig  \
    --ram_cutoff 3 --oneside_ram --min_acr 2 --min_out_conf 2 \
    >${LOGS_FILE}.job3log 2>${LOGS_FILE}.job3err


