#/bin/bash
#PBS -l walltime=72:00:00
#PBS -l nodes=1:ppn=6

# PARAMETERS

WISECONDORX_DIR="/home/leraman/tools/wisecondorX"
CORES=6

NPZ_INPUT_DIR="/home/leraman/convert.ref.output"
REF_OUTPUT_DIR="/home/leraman/newref.output"
REF_SIZES="15 50 100 200 500 1000" # space separated list (kb)
RELEASE="hg38"
GENDER="F"

# SCRIPT

for REF in ${REF_SIZES}
do
	python2 ${WISECONDORX_DIR}/wisecondorX.py newref \
	${NPZ_INPUT_DIR}/*.npz \
	${REF_OUTPUT_DIR}/reference.${RELEASE}.${REF}kb.npz \
	-binsize ${REF}000 -cpus ${CORES} -gender ${GENDER}
done