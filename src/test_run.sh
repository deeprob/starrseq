#!/bin/bash
set -ue

i_pre="/data5/deepro/starrseq/main_lib/raw_data/Input_SeqReady"
c_pre="/data5/deepro/starrseq/main_lib/raw_data/CC"
k_pre="/data5/deepro/starrseq/main_lib/raw_data/ATF2"

i_rep="A1_R1_S1 B1_R2_S2 C1_R3_S3"
c_rep="R1_STARR_Seq_S4 R2_STARR_Seq_S5 R3_STARR_Seq_S6"
k_rep="R1_STARR_Seq_S10 R2_STARR_Seq_S11 R3_STARR_Seq_S12"

i_rp="R1 R3"
c_rp="R1 R3"
k_rp="R1 R3"

ick_umi="R2 R2 R2"

i_s="001.fastq.gz"
c_s="001.fastq.gz"
k_s="001.fastq.gz"


python run_analysis.py --input_prefix $i_pre --control_prefix $c_pre --ko_prefix $k_pre --input_replicates "${i_rep}" --control_replicates "${c_rep}" --ko_replicates "${k_rep}" --input_pairs "${i_rp}" --control_pairs "${c_rp}" --ko_pairs "${k_rp}" --input_suffix $i_s --control_suffix $c_s --ko_suffix $k_s --umi "${ick_umi}" --dedup_flag --align_flag --filter_flag --bigwig_flag --peak_flag --control_flag
