PATH_BAM=$1
if [ -e ${PATH_BAM} ]; then
        BAM_SIZE=$(stat -Lc%s "$PATH_BAM")
        if (( ${BAM_SIZE} < 1 )); then
                exit 0
        fi
fi

if [ -z $2 ]; then
    TMP_PREFIX="tea"
else
    TMP_PREFIX=$2
fi

if [ -z $3 ]; then
    CORES=4
else
    CORES=$3
fi

if [ -z $4 ]; then
    MEM=50
else
    MEM=$4
fi

DIR_NAME=$(dirname $1)
BASE_NAME=$(basename $1)
PREFIX=${BASE_NAME%.bam}
TMP_DIR=${TMP_PREFIX}/${PREFIX}
FULL_PREFIX=${TMP_DIR}/${PREFIX}
mkdir -p ${FULL_PREFIX}

LOGS_DIR=${TMP_PREFIX}/logs
mkdir -p ${LOGS_DIR}
LOGS_FILE=${LOGS_DIR}/${PREFIX}

SCRIPT_ORIGINAL=02_run_tea-1file-2title-3cores.sh
SCRIPT_COPIED=${TMP_PREFIX}/02_run_tea-1file-2title-3cores.sh

cp ${SCRIPT_ORIGINAL} ${SCRIPT_COPIED} 


QUEUE="priority -t 5-0" #"medium -t 5-0" # "priority -t 0-1" # "medium -t 5-0"

sbatch -p ${QUEUE} --mail-type=ALL,TIME_LIMIT_80 -n ${CORES} --mem=${MEM}G \
	-J ${TMP_PREFIX} -o ${LOGS_FILE}.slr -e ${LOGS_FILE}.err \
	${SCRIPT_COPIED} ${PATH_BAM} ${TMP_PREFIX} ${CORES}

