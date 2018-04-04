#!/usr/bin/env bash
#PBS -l walltime=383:59:59
#PBS -l nodes=1:ppn=1

# PARAMETERS

WISECONDORX_DIR="path/to/wisecondorX" # wisecondorX clone

NPZ_FILES="path/to/samples.txt" # cases that will be tested versus the given reference
# Example of the structure of this file:
# ID_1 path/to/convert.npz/ID_1.npz
# ID_2 path/to/convert.npz/ID_2.npz
# ...
REF="path/to/newref.npz/reference.hg38.F.50kb.npz" # the cases in the NPZ_FILES document are thus expected to be female
OUTPUT_DIR="path/to/predict.output" # existing output folder


# OPTIONAL PARAMETERS

OPT="-plot -txt -json" # all output formats

# SCRIPT

while read LINE; do

    SAMPLE=$(echo $LINE | awk -F ' ' '{print $1}')
    NPZ=$(echo $LINE | awk -F ' ' '{print $2}')

    python2 ${WISECONDORX_DIR}/wisecondorX.py predict ${NPZ} ${REF} ${OUTPUT_DIR}/${SAMPLE} ${OPT}

done <${NPZ_FILES}