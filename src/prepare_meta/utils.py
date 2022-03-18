import os
import json
import numpy as np
import pandas as pd

#### This script contains all the functions that were used to create the initial data required for calling peaks from raw reads #### 

#### GLOBALS ####
CURRENT_DIR_PATH = os.path.dirname(os.path.abspath(__file__))
# the meta data dict with STARRSeq library information was created manually
# Add reference genome info and also peak caller info to the meta dict
META_DATA_DICT = {

    "input":{
        "prefix":"/data5/deepro/starrseq/lib_03162022/raw_data/IN/Input_SeqReady", 
        "replicates":"A1_R1_S1 B1_R2_S2 C1_R3_S3", 
        "read_pairs":"R1 R3", 
        "umi":"R2", 
        "suffix":"001.fastq.gz",
        "shortform":"IN"
        },
    "control":{
        "prefix":"/data5/deepro/starrseq/lib_03162022/raw_data/CC/CC", 
        "replicates":"R1_STARR_Seq_S4 R2_STARR_Seq_S5 R3_STARR_Seq_S6", 
        "read_pairs":"R1 R3", 
        "umi":"R2", 
        "suffix":"001.fastq.gz",
        "shortform":"CC"
        },
    "16p12.1":{
        "prefix":"/data5/deepro/starrseq/lib_03162022/raw_data/16P12_1/16P12_1", 
        "replicates":"R1_STARR_Seq_S25 R2_STARR_Seq_S26 R3_STARR_Seq_S27", 
        "read_pairs":"R1 R3", 
        "umi":"R2", 
        "suffix":"001.fastq.gz",
        "shortform":"16P12_1"
        },
    "atf2":{
        "prefix":"/data5/deepro/starrseq/lib_03162022/raw_data/ATF2/ATF2", 
        "replicates":"R1_STARR_Seq_S10 R2_STARR_Seq_S11 R3_STARR_Seq_S12", 
        "read_pairs":"R1 R3", 
        "umi":"R2", 
        "suffix":"001.fastq.gz",
        "shortform":"ATF2"
        },
    "ctcf":{
        "prefix":"/data5/deepro/starrseq/lib_03162022/raw_data/CTCF/CTCF", 
        "replicates":"R1_STARR_Seq_S13 R2_STARR_Seq_S14 R3_STARR_Seq_S15", 
        "read_pairs":"R1 R3", 
        "umi":"R2", 
        "suffix":"001.fastq.gz",
        "shortform":"CTCF"
        },
    "foxa1":{
        "prefix":"/data5/deepro/starrseq/lib_03162022/raw_data/FOXA1/FOXA1", 
        "replicates":"R1_STARR_Seq_S7 R2_STARR_Seq_S8 R3_STARR_Seq_S9", 
        "read_pairs":"R1 R3", 
        "umi":"R2", 
        "suffix":"001.fastq.gz",
        "shortform":"FOXA1"
        },
    "lef1":{
        "prefix":"/data5/deepro/starrseq/lib_03162022/raw_data/LEF1/LEF1", 
        "replicates":"R1_STARR_Seq_S19 R2_STARR_Seq_S20 R3_STARR_Seq_S21", 
        "read_pairs":"R1 R3", 
        "umi":"R2", 
        "suffix":"001.fastq.gz",
        "shortform":"LEF1"
        },
    "scrt1":{
        "prefix":"/data5/deepro/starrseq/lib_03162022/raw_data/SCRT1/SCRT1", 
        "replicates":"R1_STARR_Seq_S16 R2_STARR_Seq_S17 R3_STARR_Seq_S18", 
        "read_pairs":"R1 R3", 
        "umi":"R2", 
        "suffix":"001.fastq.gz",
        "shortform":"SCRT1"
        },
    "tcf7l2":{
        "prefix":"/data5/deepro/starrseq/lib_03162022/raw_data/TCF7L2/TCF7L2", 
        "replicates":"R1_STARR_Seq_S22 R2_STARR_Seq_S23 R3_STARR_Seq_S24", 
        "read_pairs":"R1 R3", 
        "umi":"R2", 
        "suffix":"001.fastq.gz",
        "shortform":"TCF7L2"
        },
    "genome":{
        "ref_fasta":"/data5/deepro/genomes/hg38/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta",
        "ref_twobit":"/data5/deepro/genomes/hg38/hg38.2bit"
        },
    "roi":{
        "roi_sorted":"/data5/deepro/starrseq/computational_pipeline/data/roi/master.sorted.bed"
        },
    "starrpeaker": {
        "data_dir":"/data5/deepro/starrseq/starrpeaker_data/"
    },
    "cradle": {
        "data_dir":"/data5/deepro/starrseq/cradle_data/"
    },

}

#######################################
# master list preprocessing functions #
#######################################

# xlsx to bed with TF annotations function
def get_tfs(df_row):
    return list(set(df_row[6:].loc[~df_row[6:].isna()]))

def excel_to_bed_with_annot(in_excel, out_bed):
    df_excel = pd.read_excel(in_excel)
    # get the tfs per chromosomal location as a list
    df_excel["tfs"] = df_excel.apply(get_tfs, axis=1)
    # only keep chrm, start, end, size, tfs columns
    df_excel = df_excel.iloc[:, [0,1,2,4,-1]]
    # get all the tfs that are annotation along all locations
    all_tfs = list(sorted(set(sum(df_excel.tfs.values.flatten(), []))))
    # create a dictionary that maps sorted tfs to a number
    tf_dict = dict(zip(all_tfs, range(len(all_tfs))))
    # create an array for one hot encoding the tfs per location
    tf_array = np.zeros((len(df_excel), len(all_tfs)))
    for idx, tf_list in enumerate(df_excel.tfs):
        for tf in tf_list:
            tf_array[idx, tf_dict[tf]] = 1
    # concatenate the array with the original df
    df_tf = pd.concat((df_excel.iloc[:, :4], pd.DataFrame(tf_array, columns=all_tfs)), axis=1)
    # save the array to disk
    df_tf.to_csv(out_bed, sep="\t", index=False)
    return 


# TODO: bed sorted by chromosomes

# TODO: fa file generated for chromosomal locations - this file is not required right now

# TODO: chromosomal locations with sequence in bed format - this file is also not required


#####################
# metadata creation #
#####################

def create_metadata(metadata_filename):
    """
    Create a metadata dict in json format that contains info about all the 
    STARRSeq libraries and store them in the REPO_HOME/data/lib_meta
    """

    with open(metadata_filename, "w") as f:  
        json.dump(META_DATA_DICT, f, indent=4)
    return
