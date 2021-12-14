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
conda activate cradle

# get all arguments
while getopts i:r: flag
do
    case "${flag}" in
        i) INPUT_PREFIX=${OPTARG};;
        r) BIOL_REPS=(${OPTARG});;
    esac
done

# example command: 
# bash bamtobw.sh /data5/deepro/starrseq/main_lib/aligned_reads/Input_SeqReady_C1_R3_S3_filtered.bam /data5/deepro/starrseq/main_lib/aligned_reads/Input_SeqReady_C1_R3_S3_filtered.bw

for rep in "${BIOL_REPS[@]}" # ${AR[0]}
do
    bamfile=${INPUT_PREFIX}_${rep}_filtered.bam
    outfile=${INPUT_PREFIX}_${rep}_filtered.bw

    samtools index $bamfile
    bamCoverage -b $bamfile -o $outfile -of bigwig -p 64
done
