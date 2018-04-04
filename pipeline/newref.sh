#!/usr/bin/env bash
#PBS -l walltime=383:59:59
#PBS -l nodes=1:ppn=6

# PARAMETERS

CORES=6
WISECONDORX_DIR="path/to/wisecondorX" # wisecondorX clone

INPUT_DIR="path/to/convert.npz" # existing (non-empty) folder, containing reference .npz files
OUTPUT_DIR="path/to/newref.npz" # existing (empty) folder
REF_SIZES="15 50 100 200 500 1000" # space separated list (kb)
RELEASE="hg38" # reference used to create bam files (solely used for reference filename)
GENDER="F" # gender of the of the cases at NPZ_INPUT_DIR

# SCRIPT

for REF in ${REF_SIZES}
do
	python2 ${WISECONDORX_DIR}/wisecondorX.py newref \
	${INPUT_DIR}/*.npz \
	${OUTPUT_DIR}/reference.${RELEASE}.${GENDER}.${REF}kb.npz \
	-binsize ${REF}000 -cpus ${CORES} -gender ${GENDER}
done