#!/bin/bash
set -ue

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/afs/bx.psu.edu/user/d/dzb5732/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/afs/bx.psu.edu/user/d/dzb5732/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/afs/bx.psu.edu/user/d/dzb5732/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/afs/bx.psu.edu/user/d/dzb5732/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup

# activate conda environment
conda activate starrpeaker_db

## Get file paths to required arguments ##
# Prefix to store peaks
prefix=$1

# starrseq input file
input=$2

#starrseq output file
output=$3

# chromosome sizes file
chromsize="/data5/deepro/starrseq/starrpeaker_data/GRCh38.chrom.sizes.simple.sorted" 

# blacklist regions file
blacklist="/data5/deepro/starrseq/starrpeaker_data/ENCODE_blacklist_GRCh38_ENCFF419RSJ_merged.bed"

# covariate files
cov1="/data5/deepro/starrseq/covariates/STARRPeaker_cov_GRCh38_gem-mappability-100mer.bw"
cov2="/data5/deepro/starrseq/covariates/STARRPeaker_cov_GRCh38_ucsc-gc-5bp.bw"
cov3="/data5/deepro/starrseq/covariates/STARRPeaker_cov_GRCh38_linearfold-folding-energy-100bp.bw"

# starrpeaker peak calling command
starrpeaker --prefix $prefix --chromsize $chromsize --blacklist $blacklist --cov $cov1 $cov2 $cov3 -i $input -o $output
