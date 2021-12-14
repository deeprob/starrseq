#!/bin/bash
set -ue

# This script aligns all the fq reads of a library to the human reference genome
# assumptions: 
# 0) environment installed and activated with proper software
# 1) ref genome is downloaded
# 2) ref genome is indexed using bwa index and the indexed files are in the same dir as the ref genome
# 3) Naming convention of the files to be aligned: {INPUT_PREFIX}_{REPLICATE_INDEX}_{PAIR_INDEX}_{SUFFIX}

# example command
# $ bash align_fastq_to_ref.sh -g "/data5/deepro/genomes/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta" -i "/data5/Moi/MiSeq_Test_fullPool_1/Input_SeqReady_A1_R1" -r "S1" -p "R1 R3" -s "001" -o "/data5/deepro/starrseq/miseq_test/aligned_reads/input"

# get all arguments
while getopts g:i:r:p:s:o: flag
do
    case "${flag}" in
        g) REF_GENOME=${OPTARG};;
        i) INPUT_PREFIX=${OPTARG};;
        r) BIOL_REPS=(${OPTARG});;
        p) PAIRS=(${OPTARG});;
        s) SUFFIX=${OPTARG};;        
        o) OUTPUT_PREFIX=${OPTARG};;

    esac
done

# align reads throughout the library

for rep in "${BIOL_REPS[@]}" # ${AR[0]}
do
    # align reads
    bwa mem -t 64 -v 2 ${REF_GENOME} ${INPUT_PREFIX}_${rep}_${PAIRS[0]}_${SUFFIX} ${INPUT_PREFIX}_${rep}_${PAIRS[1]}_${SUFFIX} > ${OUTPUT_PREFIX}_${rep}_unfiltered.sam
    # convert sam to bam
    samtools sort -@ 24 -o ${OUTPUT_PREFIX}_${rep}_unfiltered.bam ${OUTPUT_PREFIX}_${rep}_unfiltered.sam 
    # delete the sam file, takes up too much space
    rm ${OUTPUT_PREFIX}_${rep}_unfiltered.sam
done
