#!/bin/bash
set -ue

roi_file=$1
input_library_filtered_bam=$2
input_library_coverage_bed=$3
depth=$4

if [ $depth == true ]
then
    bedtools coverage -a $roi_file -b $input_library_filtered_bam -d > $input_library_coverage_bed
else
    bedtools coverage -a $roi_file -b $input_library_filtered_bam > $input_library_coverage_bed
fi
