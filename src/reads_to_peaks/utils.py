import os
import subprocess
import json
from argparse import Namespace

# internal imports : since I am running it from the top level
import starrdust.starrdust as sd


#### This script contains all the functions that were used to create the intermediate data required for calling peaks from raw reads #### 

#### GLOBALS ####
CURRENT_DIR_PATH = os.path.dirname(os.path.abspath(__file__))

#########################
# args list for running #
#########################

def create_args(meta_file, ko_name,
    input_flag=True, output_flag=True, 
    dedup_flag=True, align_flag=True, filter_flag=True, 
    starrpeaker_peak_flag=True, cradle_peak_flag=True, macs2_peak_flag=True):

    with open(meta_file, "r") as f: 
        meta_dict = json.load(f)

    args = Namespace(
        # from metadata and globals
        input_library_prefix = meta_dict["input"]["prefix"],
        output_library_prefix = meta_dict[ko_name]["prefix"],
        input_library_reps = meta_dict["input"]["replicates"],
        output_library_reps = meta_dict[ko_name]["replicates"],
        input_library_pair = meta_dict["input"]["read_pairs"],
        output_library_pair= meta_dict[ko_name]["read_pairs"],
        input_library_umi = meta_dict["input"]["umi"],
        output_library_umi = meta_dict[ko_name]["umi"],
        input_library_suffix = meta_dict["input"]["suffix"],
        output_library_suffix = meta_dict[ko_name]["suffix"],
        input_library_short = meta_dict["input"]["shortform"],
        output_library_short = meta_dict[ko_name]["shortform"],
        reference_genome = meta_dict["genome"]["ref_fasta"],
        reference_genome_twobit = meta_dict["genome"]["ref_twobit"],
        roi_file = meta_dict["roi"]["roi_sorted"],
        starrpeaker_data_dir = meta_dict["starrpeaker"]["data_dir"],
        cradle_data_dir = meta_dict["cradle"]["data_dir"],
        # TODO: store short form info as well to create path later on!

        # from global flags
        input_flag = input_flag,
        output_flag = output_flag,
        dedup_flag = dedup_flag,
        align_flag = align_flag,
        filter_flag = filter_flag,
        starrpeaker_peak_flag = starrpeaker_peak_flag,
        cradle_peak_flag = cradle_peak_flag,
        macs2_peak_flag = macs2_peak_flag,

    )
    return args


#########################
# deduplication helpers #
#########################

# starrdust helpers
def get_deduped_files_helper(read1, read2, umi, mapq_r1=30, mapq_r2=30, min_length_r1=100, min_length_r2=100):
    sd.starrdust(read1, read2, umi, mapq_r1, mapq_r2, min_length_r1, min_length_r2, None)
    return

def get_deduped_files(lib_pre, lib_rep, lib_pairs, lib_umi_idx, lib_suff):
    lib_rep_list = lib_rep.split()
    lib_pair_list = lib_pairs.split()

    for rep in lib_rep_list:
        lib_read1_path = "_".join([lib_pre, rep, lib_pair_list[0], lib_suff])
        lib_read2_path = "_".join([lib_pre, rep, lib_pair_list[1], lib_suff])
        lib_umi_path = "_".join([lib_pre, rep, lib_umi_idx, lib_suff])

        get_deduped_files_helper(lib_read1_path, lib_read2_path, lib_umi_path)
    return


#####################
# alignment helpers #
#####################

def align_reads_helper(reference_genome, library_prefix, library_replicates, 
                       library_read_pairs, library_suffix, library_aligned_prefix, 
                       cores=64):
    """
    Align reads using bwa - bash script under the ./shell_scripts dir
    """
    os.makedirs(os.path.dirname(library_aligned_prefix), exist_ok=True)
    cmd = [
        "bash", f"{CURRENT_DIR_PATH}/shell_scripts/1_align_fastq_to_ref.sh", 
        "-g", f"{reference_genome}", "-i", f"{library_prefix}", 
        "-r", f"{library_replicates}", "-p", f"{library_read_pairs}", 
        "-s", f"{library_suffix}", "-o", f"{library_aligned_prefix}",
         "-t", f"{cores}"]

    subprocess.run(cmd)
    return


##################
# filter helpers #
##################

def filter_reads_helper(library_prefix, library_replicates, library_filtered_prefix, 
                        roi_file, dedup_flag):
    """
    Filter reads using samtools - bash script under the ./shell_scripts dir
    TODO: Think about transferring this to pysam since it does not need a new environment and has a python tool
    """
    os.makedirs(os.path.dirname(library_filtered_prefix), exist_ok=True)
    cmd = [
        "bash", f"{CURRENT_DIR_PATH}/shell_scripts/2_filter_reads.sh", 
        "-i", f"{library_prefix}", "-r", f"{library_replicates}", 
        "-o", f"{library_filtered_prefix}", "-f", f"{roi_file}", 
        "-d", f"{dedup_flag}"]

    subprocess.run(cmd)
    return


########################
# peak calling helpers #
########################

# starrpeaker peak calling functions
def call_starrpeaker_peaks_helper(peaks_prefix, input_filtered_bam, output_filtered_bam, starrpeaker_data_dir):

    cmd = ["bash", f"{CURRENT_DIR_PATH}/shell_scripts/3a_call_peaks_starrpeaker.sh", 
            f"{peaks_prefix}", 
            f"{input_filtered_bam}", f"{output_filtered_bam}", 
            f"{starrpeaker_data_dir}"]
    subprocess.run(cmd)
    return

def call_starrpeaker_peaks(input_library_filtered_prefix, 
                        output_library_filtered_prefix,
                        starrpeaker_data_dir,
                        output_flag=False):
    """Call peaks for output library using starrpeaker"""

    # create the output peaks path where it will be stored
    output_peaks_prefix = os.path.join(os.path.dirname(output_library_filtered_prefix.replace("filtered_libraries", "results/peaks")), "starrpeaker")
    # create the output peaks path where it will be stored
    os.makedirs(output_peaks_prefix, exist_ok=True)
    # prepare the input files, starrpeaker requires merged bam
    input_library_filtered_bam =  input_library_filtered_prefix + ".bam"
    output_library_filtered_bam =  output_library_filtered_prefix + ".bam"

    if output_flag:
        output_peaks_prefix = os.path.join(output_peaks_prefix, "peaks")
        call_starrpeaker_peaks_helper(output_peaks_prefix, input_library_filtered_bam, output_library_filtered_bam, starrpeaker_data_dir)

    return

# cradle peak calling functions
def create_bigwig_helper(library_prefix, library_replicates):
    cmd = ["bash", f"{CURRENT_DIR_PATH}/shell_scripts/3b_0_bamtobw.sh", "-i", f"{library_prefix}", "-r", f"{library_replicates}"]
    subprocess.run(cmd)
    return

def create_bigwig(input_library_prefix, input_library_replicates,
                  output_library_prefix, output_library_replicates,
                  input_flag=False, output_flag=False):
    """
    Create bigwig files from the aligned and filtered bamfiles for visualization and peak calling using cradle
    """
    if input_flag:
        create_bigwig_helper(input_library_prefix, input_library_replicates)
    if output_flag:
        create_bigwig_helper(output_library_prefix, output_library_replicates)

    return

def call_cradle_peaks_helper(input_rep_dir, output_rep_dir, 
                             input_rep_paths, output_rep_paths, 
                             peaks_prefix, 
                             roi_file, reference_genome_twobit,
                             cradle_data_dir):
    
    cmd = ["bash", f"{CURRENT_DIR_PATH}//shell_scripts/3b_1_call_peaks_cradle.sh"]
    cmd += [input_rep_dir, output_rep_dir, peaks_prefix]
    cmd += input_rep_paths
    cmd += output_rep_paths
    cmd += [roi_file, reference_genome_twobit, cradle_data_dir]
    subprocess.run(cmd)
    return

def call_cradle_peaks(input_library_prefix, input_library_replicates,
                      output_library_prefix, output_library_replicates,
                      reference_genome_twobit, roi_file, 
                      cradle_data_dir,   
                      output_flag=False):

    """
    Call peaks for output library using cradle
    """

    output_peaks_prefix = os.path.join(os.path.dirname(output_library_prefix.replace("aligned_reads", "results/peaks")), "cradle")
    # create the output peaks path where it will be stored
    os.makedirs(output_peaks_prefix, exist_ok=True)
    # get the paths where the replicates are stored from the library prefixes
    input_rep_dir = os.path.dirname(input_library_prefix)
    input_rep_paths = list(map(lambda x: os.path.basename("_".join([input_library_prefix, x, "filtered"])) , input_library_replicates.split(" ")))
    output_rep_dir = os.path.dirname(output_library_prefix)
    output_rep_paths = list(map(lambda x: os.path.basename("_".join([output_library_prefix, x, "filtered"])) , output_library_replicates.split(" ")))

    if output_flag:
        call_cradle_peaks_helper(input_rep_dir, output_rep_dir, 
                             input_rep_paths, output_rep_paths, 
                             output_peaks_prefix, roi_file, reference_genome_twobit,
                             cradle_data_dir)

    return

# macs2 peak calling functions
def call_macs2_peaks_helper(peaks_prefix, input_filtered_bam, output_filtered_bam):

    cmd = ["bash", f"{CURRENT_DIR_PATH}/shell_scripts/3c_call_peaks_macs2.sh", 
            f"{peaks_prefix}", 
            f"{input_filtered_bam}", f"{output_filtered_bam}"]
    subprocess.run(cmd)
    return

def call_macs2_peaks(input_library_filtered_prefix, 
                    output_library_filtered_prefix,
                    output_flag=False):
    """
    Call peaks for output library using macs2
    """

    # create the output peaks path where it will be stored
    output_peaks_prefix = os.path.join(os.path.dirname(output_library_filtered_prefix.replace("filtered_libraries", "results/peaks")), "macs2")
    # create the output peaks path where it will be stored
    os.makedirs(output_peaks_prefix, exist_ok=True)
    # prepare the input files, macs2 requires merged bam
    input_library_filtered_bam =  input_library_filtered_prefix + ".bam"
    output_library_filtered_bam =  output_library_filtered_prefix + ".bam"

    if output_flag:
        call_macs2_peaks_helper(output_peaks_prefix, input_library_filtered_bam, output_library_filtered_bam)

    return
