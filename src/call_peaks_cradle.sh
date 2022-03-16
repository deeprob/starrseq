#!/bin/bash
set -ue

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/data5/deepro/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/data5/deepro/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/data5/deepro/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/data5/deepro/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

# activate conda environment
conda activate cradle

## Get file paths to required arguments ##
# Dir to store bias corrected bigwigs
outprefix_corrected=$1
# Dir to store peaks
outprefix_peaks=$2

# starrseq input file prefix
input1=$3
input2=$4
input3=$5

#starrseq output files
output1=$6
output2=$7
output3=$8

# roi_file
roi_file=$9
# genome_file
genome_file=${10}

# echo $outprefix_corrected
# echo $outprefix_peaks
# echo $input1
# echo $input2
# echo $input3
# echo $output1
# echo $output2
# echo $output3
# echo $roi_file
# echo $genome_file

# blacklist regions file
blacklist="/data5/deepro/starrseq/starrpeaker_data/ENCODE_blacklist_GRCh38_ENCFF419RSJ_merged.bed"

# mappability bias file
mapfile="/data5/deepro/starrseq/cradle_data/mappability/hg38_mappability_100mer.bw" 
# gquad bias file
gquadfile="/data5/deepro/starrseq/cradle_data/gquadruplex/GSE63874_Na_K_PDS_plus_hits_intersect_hg38_uniq_K.bw"

# covariate directory stored file
cov_dir="/data5/deepro/starrseq/cradle_data/hg38_fragLen500_kmer100"

# cradle correctBias command
cradle correctBias -ctrlbw ${input1}.bw ${input2}.bw ${input3}.bw -expbw ${output1}.bw ${output2}.bw ${output3}.bw -l 500 -r ${roi_file} -biasType shear pcr map gquad -genome ${genome_file} -bl $blacklist -p 64 -o $outprefix_corrected -kmer 100 -mapFile $mapfile -gquadFile $gquadfile
# cradle correctBias_stored command
# cradle correctBias_stored -ctrlbw ${input1}.bw ${input2}.bw ${input3}.bw -expbw ${output1}.bw ${output2}.bw ${output3}.bw -r ${roi_file} -biasType shear pcr map gquad -genome ${genome_file} -bl $blacklist -p 64 -o $outprefix_corrected -covariDir $cov_dir

# # cradle peakcall command
cradle callPeak -ctrlbw ${input1}_corrected.bw ${input2}_corrected.bw ${input3}_corrected.bw -expbw ${output1}_corrected.bw ${output2}_corrected.bw ${output3}_corrected.bw -r $roi_file -fdr 0.05 -o $outprefix_peaks -p 64


## command 
# bash call_peaks_cradle.sh /data5/deepro/starrseq/main_lib/aligned_reads /data5/deepro/starrseq/main_lib/results/peaks/CC '/data5/deepro/starrseq/main_lib/aligned_reads/Input_SeqReady_A1_R1_S1_filtered' '/data5/deepro/starrseq/main_lib/aligned_reads/Input_SeqReady_B1_R2_S2_filtered' '/data5/deepro/starrseq/main_lib/aligned_reads/Input_SeqReady_C1_R3_S3_filtered' '/data5/deepro/starrseq/main_lib/aligned_reads/CC_R1_STARR_Seq_S4_filtered' '/data5/deepro/starrseq/main_lib/aligned_reads/CC_R2_STARR_Seq_S5_filtered' '/data5/deepro/starrseq/main_lib/aligned_reads/CC_R3_STARR_Seq_S6_filtered' '/afs/bx.psu.edu/user/d/dzb5732/work/girirajan_lab/starrseq/data/master.sorted.bed' '/data5/deepro/genomes/hg38.2bit'
