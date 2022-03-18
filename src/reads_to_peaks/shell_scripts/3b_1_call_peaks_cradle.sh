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
# dir that points to the starrseq input bigwig files
indir=$1
# dir that points to the starrseq output bigwig files
outdir=$2
# dir that points to the output files where bias corrected bw and peaks will be stored
peakdir=$3

# starrseq input file basenames
input1=$4
input2=$5
input3=$6

#starrseq output files basenames
output1=$7
output2=$8
output3=$9

# roi_file
roi_file=${10}
# genome_file
genome_file=${11}
# cradle data directory that stores covariates
cradle_dir=${12}

# blacklist regions file
blacklist="${cradle_dir}ENCODE_blacklist_GRCh38_ENCFF419RSJ_merged.bed"

# mappability bias file
mapfile="${cradle_dir}mappability/hg38_mappability_100mer.bw" 
# gquad bias file
gquadfile="${cradle_dir}gquadruplex/GSE63874_Na_K_PDS_plus_hits_intersect_hg38_uniq_K.bw"

# covariate directory stored file
cov_dir="${cradle_dir}hg38_fragLen500_kmer100"

# cradle correctBias command
cradle correctBias -ctrlbw ${indir}/${input1}.bw ${indir}/${input2}.bw ${indir}/${input3}.bw -expbw ${outdir}/${output1}.bw ${outdir}/${output2}.bw ${outdir}/${output3}.bw -l 500 -r ${roi_file} -biasType shear pcr map gquad -genome ${genome_file} -bl ${blacklist} -p 64 -o ${peakdir} -kmer 100 -mapFile $mapfile -gquadFile $gquadfile
# cradle correctBias_stored command
# cradle correctBias_stored -ctrlbw ${input1}.bw ${input2}.bw ${input3}.bw -expbw ${output1}.bw ${output2}.bw ${output3}.bw -r ${roi_file} -biasType shear pcr map gquad -genome ${genome_file} -bl $blacklist -p 64 -o $outprefix_corrected -covariDir $cov_dir

# # cradle peakcall command
cradle callPeak -ctrlbw ${peakdir}/${input1}_corrected.bw ${peakdir}/${input2}_corrected.bw ${peakdir}/${input3}_corrected.bw -expbw ${peakdir}/${output1}_corrected.bw ${peakdir}/${output2}_corrected.bw ${peakdir}/${output3}_corrected.bw -r $roi_file -fdr 0.05 -o ${peakdir} -p 64


## command 
# bash call_peaks_cradle.sh /data5/deepro/starrseq/main_lib/aligned_reads/IN /data5/deepro/starrseq/main_lib/aligned_reads/CC /data5/deepro/starrseq/main_lib/results/peaks/CC 'Input_SeqReady_A1_R1_S1_filtered' 'Input_SeqReady_B1_R2_S2_filtered' 'Input_SeqReady_C1_R3_S3_filtered' 'CC_R1_STARR_Seq_S4_filtered' 'CC_R2_STARR_Seq_S5_filtered' 'CC_R3_STARR_Seq_S6_filtered' '/data5/deepro/starrseq/computational_pipeline/data/roi/master.sorted.bed' '/data5/deepro/genomes/hg38.2bit'
