#!/bin/bash
set -ue

i_pre="/data5/deepro/starrseq/nextseq_feb_2021/raw_data/SCP1CompleteLibrary"
c_pre="/data5/deepro/starrseq/nextseq_feb_2021/raw_data/SCP1_293_CC"
k_pre="/data5/deepro/starrseq/nextseq_feb_2021/raw_data/CP1_293_ATF2"

i_rep="1_S4 2_S5 3_S6"
c_rep="1_S21 2_S22 3_S23"
k_rep="1_S24 2_S25 3_S26"

i_rp="R1 R2"
c_rp="R1 R2"
k_rp="R1 R2"

ick_umi=""

i_s="001.fastq.gz"
c_s="001.fastq.gz"
k_s="001.fastq.gz"


python run_analysis.py --input_prefix $i_pre --control_prefix $c_pre --ko_prefix $k_pre --input_replicates "${i_rep}" --control_replicates "${c_rep}" --ko_replicates "${k_rep}" --input_pairs "${i_rp}" --control_pairs "${c_rp}" --ko_pairs "${k_rp}" --input_suffix $i_s --control_suffix $c_s --ko_suffix $k_s --umi "${ick_umi}" --dedup_flag --align_flag 
