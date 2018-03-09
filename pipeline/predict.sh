#/bin/bash
#PBS -l walltime=72:00:00
#PBS -l nodes=1:ppn=1

# PARAMETERS

WISECONDORX_DIR="/home/leraman/tools/wisecondorX"

NPZs_FILE="/Users/leraman/samples.txt"
# Example of the structure of this file:
# D4984165 x/x/convert.case.output/D4984165.npz
# D4984166 x/x/convert.case.output/D4984166.npz
# ...
REF="/Users/leraman/newref.output/reference.hg38.50kb.npz"
OUTPUT_DIR="/Users/leraman/predict.output"


# OPTIONAL PARAMETERS

OPT="-plot"

# SCRIPT

while read LINE; do

    SAMPLE=$(echo $LINE | awk -F ' ' '{print $1}')
    NPZ=$(echo $LINE | awk -F ' ' '{print $2}')

    python2 ${WISECONDORX_DIR}/wisecondorX.py predict ${NPZ} ${REF} ${OUTPUT_DIR}/${SAMPLE} ${OPT}

done <${NPZs_FILE}