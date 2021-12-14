#!/bin/bash
set -ue

region_file=$1
reference_genome=$2
background_region=$3
output_file=$4

# TODO: add background bed file
findMotifsGenome.pl $region_file $reference_genome $output_file -bg $background_region -p 64