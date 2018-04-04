#!/usr/bin/env bash
#PBS -l walltime=383:59:59
#PBS -l nodes=1:ppn=1

# PARAMETERS

WISECONDORX_DIR="path/to/wisecondorX" # wisecondorX clone

BAM_FILES="path/to/bam_files.txt" # reference or test cases
# Example of the structure of this file:
# ID_1 path/to/ID_1.bam
# ID_2 path/to/ID_2.bam
OUTPUT_DIR="path/to/convert.npz"

# SCRIPT

while read LINE; do

    SAMPLE=$(echo $LINE | awk -F ' ' '{print $1}')
    BAM=$(echo $LINE | awk -F ' ' '{print $2}')

    #Create bins @5kb
    echo "creating 5kb bins for sample ${SAMPLE}"
    python2 ${WISECONDORX_DIR}/wisecondorX.py convert ${BAM} ${OUTPUT_DIR}/${SAMPLE}.npz

done < ${BAM_FILES}
