#!/bin/bash
set -ue

major_file=$1
minor_file=$2
output_file=$3
common_flag=$4

if [ $common_flag == true ]
then
    # get set intersecting elements of the major file with the minor file 
    bedtools intersect -a $major_file -b $minor_file > $output_file
else
    # get the non intersecting elements of the major file with the minor file
    bedtools intersect -v -a $major_file -b $minor_file > $output_file
fi
