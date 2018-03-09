#/bin/bash
#PBS -l walltime=383:59:59
#PBS -l nodes=1:ppn=4

# PARAMETERS

WISECONDORX_DIR="/home/leraman/tools/wisecondorX"
CORES=4

FASTQ_FILES="/home/leraman/fastq.txt"
# Example of the structure of this file:
# D4984165 ././D4984165_S05_xxx_002.fastq.gz,././D4984165_S05_xxx_003.fastq.gz,././D4984165_S05_xxx_004.fastq.gz
# D4984166 ././D4984165_S06_xxx_002.fastq.gz,././D4984166_S05_xxx_003.fastq.gz,././D4984166_S05_xxx_004.fastq.gz
BOWTIE2_INDEX="/Shared/references/Hsapiens/hg38/hg38/hg38full/bowtie2_index/hg38full"
NPZ_OUTPUT_DIR="/home/leraman/convert.ref.output"

# SCRIPT

cd ${NPZ_OUTPUT_DIR}

while IFS='' read -r LINE || [[ -n "$line" ]]; do

    SAMPLE=$(echo $LINE | awk -F ' ' '{print $1}')
    FASTQS=$(echo $LINE | awk -F ' ' '{print $2}')

    BOWTIE_INPUT=""
    for FASTQ in $(echo ${FASTQS} | sed "s/,/ /g") ; do
        BOWTIE_INPUT="${BOWTIE_INPUT},${FASTQ}"
    done
    BOWTIE_INPUT=$(echo ${BOWTIE_INPUT:1})

    #Map fastqs
    echo "mapping sample ${SAMPLE}"

    bowtie2 --local \
      -p ${CORES} \
      --fast-local \
      -x ${BOWTIE2_INDEX} \
      -U ${BOWTIE_INPUT} \
    | bamsormadup \
      inputformat=sam \
      threads=${CORES} \
      SO=coordinate \
      outputformat=bam \
      indexfilename=${SAMPLE}.bam.bai > ${SAMPLE}.bam

    #Create bins @5kb
    echo "creating 5kb bins for sample ${SAMPLE}"
    python2 ${WISECONDORX_DIR}/wisecondorX.py convert ${SAMPLE}.bam ${SAMPLE}.npz
    rm ${SAMPLE}.bam.bai
    rm ${SAMPLE}.bam

done < ${FASTQ_FILES}
