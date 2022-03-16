# Import modules
import os
import subprocess
from pprint import pprint

# internal imports
from starrdust.starrdust import starrdust

CURRENT_DIR_PATH = os.path.dirname(os.path.abspath(__file__))

# TODO: make sure that all files exists using the following function
def file_exists(filepath):
    return os.path.exists(filepath)

# starrdust helpers
def get_deduped_files_helper(read1, read2, umi, mapq_r1=30, mapq_r2=30, min_length_r1=100, min_length_r2=100):
    starrdust(read1, read2, umi, mapq_r1, mapq_r2, min_length_r1, min_length_r2, None)
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

# remove duplicates with starrdust
def remove_dups(input_library_prefix, input_library_replicates, input_library_read_pairs, input_library_umi_index, input_library_suffix,
                control_library_prefix, control_library_replicates, control_library_read_pairs, control_library_umi_index, control_library_suffix,
                ko_library_prefix, ko_library_replicates, ko_library_read_pairs, ko_library_umi_index, ko_library_suffix, 
                input_flag=False, control_flag=False):
    """Remove duplicates from all three libraries using starrdust"""

    if input_flag:
        # input library deduplication
        get_deduped_files(input_library_prefix, input_library_replicates, input_library_read_pairs, input_library_umi_index, input_library_suffix)
    if control_flag:
        # control library deduplication
        get_deduped_files(control_library_prefix, control_library_replicates, control_library_read_pairs, control_library_umi_index, control_library_suffix)
    # ko library deduplication
    get_deduped_files(ko_library_prefix, ko_library_replicates, ko_library_read_pairs, ko_library_umi_index, ko_library_suffix)
    return

# align reads to the reference genome
def align_reads(input_library_prefix, input_library_replicates, input_library_read_pairs, input_library_suffix, input_library_aligned_prefix,
                control_library_prefix, control_library_replicates, control_library_read_pairs, control_library_suffix, control_library_aligned_prefix,
                ko_library_prefix, ko_library_replicates, ko_library_read_pairs, ko_library_suffix, ko_library_aligned_prefix, 
                reference_genome, input_flag=False, control_flag=False):
    """Align reads to the reference genome"""

    if input_flag:
        # align input reads to the reference genome
        subprocess.call(["bash", f"{CURRENT_DIR_PATH}/align_fastq_to_ref.sh", "-g", f"{reference_genome}", "-i", f"{input_library_prefix}", "-r", f"{input_library_replicates}", "-p", f"{input_library_read_pairs}", "-s", f"{input_library_suffix}", "-o", f"{input_library_aligned_prefix}"])
    if control_flag:
        # align control reads to the reference genome
        subprocess.call(["bash", f"{CURRENT_DIR_PATH}/align_fastq_to_ref.sh", "-g", f"{reference_genome}",  "-i", f"{control_library_prefix}", "-r", f"{control_library_replicates}", "-p", f"{control_library_read_pairs}", "-s", f"{control_library_suffix}", "-o", f"{control_library_aligned_prefix}"])
    # align ko reads to the reference genome
    subprocess.call(["bash", f"{CURRENT_DIR_PATH}/align_fastq_to_ref.sh", "-g", f"{reference_genome}",  "-i", f"{ko_library_prefix}", "-r", f"{ko_library_replicates}", "-p", f"{ko_library_read_pairs}", "-s", f"{ko_library_suffix}", "-o", f"{ko_library_aligned_prefix}"])
    return

# filter reads
def filter_reads(input_library_aligned_prefix, input_library_replicates, input_library_filtered_prefix,
                control_library_aligned_prefix, control_library_replicates, control_library_filtered_prefix,
                ko_library_aligned_prefix, ko_library_replicates, ko_library_filtered_prefix, roi_file,
                umi=True, input_flag=False, control_flag=False):
    """filter bad reads"""
    dedup = "true"
    if umi:
        dedup = "false"

    if input_flag:
        # filter input reads samtools -F 2828 or 2852 (depending on the use of starrdust or picard), -f 2 -q 30
        subprocess.call(["bash", f"{CURRENT_DIR_PATH}/filter_reads.sh", "-i", f"{input_library_aligned_prefix}", "-r", f"{input_library_replicates}", "-o", f"{input_library_filtered_prefix}", "-f", f"{roi_file}", "-d", f"{dedup}"])
    if control_flag:
        # filter control reads
        subprocess.call(["bash", f"{CURRENT_DIR_PATH}/filter_reads.sh", "-i", f"{control_library_aligned_prefix}", "-r", f"{control_library_replicates}", "-o", f"{control_library_filtered_prefix}", "-f", f"{roi_file}", "-d", f"{dedup}"])
    # filter ko reads 
    subprocess.call(["bash", f"{CURRENT_DIR_PATH}/filter_reads.sh", "-i", f"{ko_library_aligned_prefix}", "-r", f"{ko_library_replicates}", "-o", f"{ko_library_filtered_prefix}", "-f", f"{roi_file}", "-d", f"{dedup}"])
    return 

def call_peaks_helper(input_library_filtered_prefix, control_library_filtered_prefix, ko_library_filtered_prefix, control_flag=False):
    """Call peaks for control and ko library using starrpeaker"""
    control_peaks_prefix = os.path.join(control_library_filtered_prefix.replace("filtered_libraries", "results/peaks"), "starrpeaker")
    ko_peaks_prefix = os.path.join(ko_library_filtered_prefix.replace("filtered_libraries", "results/peaks"), "starrpeaker")

    input_library_filtered_bam =  input_library_filtered_prefix + ".bam"
    control_library_filtered_bam =  control_library_filtered_prefix + ".bam"
    ko_library_filtered_bam =  ko_library_filtered_prefix + ".bam"

    if control_flag:
        os.makedirs(control_peaks_prefix, exist_ok=True)
        subprocess.call(["bash", f"{CURRENT_DIR_PATH}/call_peaks.sh", f"{control_peaks_prefix}", f"{input_library_filtered_bam}", f"{control_library_filtered_bam}"])

    os.makedirs(ko_peaks_prefix, exist_ok=True)
    subprocess.call(["bash", f"{CURRENT_DIR_PATH}/call_peaks.sh", f"{ko_peaks_prefix}", f"{input_library_filtered_bam}", f"{ko_library_filtered_bam}"])
    return

def create_bigwig(input_library_prefix, input_library_replicates,
                  control_library_prefix, control_library_replicates,
                  ko_library_prefix, ko_library_replicates,
                  input_flag=False, control_flag=False):
    """Create bigwig files from the aligned and filtered bamfiles for visualization and peak calling using cradle"""
    if input_flag:
        subprocess.call(["bash", f"{CURRENT_DIR_PATH}/bamtobw.sh", "-i", f"{input_library_prefix}", "-r", f"{input_library_replicates}"])
    if control_flag:
        subprocess.call(["bash", f"{CURRENT_DIR_PATH}/bamtobw.sh",  "-i", f"{control_library_prefix}", "-r", f"{control_library_replicates}"])

    subprocess.call(["bash", f"{CURRENT_DIR_PATH}/bamtobw.sh", "-i", f"{ko_library_prefix}", "-r", f"{ko_library_replicates}"])
    return

# TODO: call peaks helper cradle    
def call_peaks_helper_cradle(input_library_prefix, input_library_replicates,
                             control_library_prefix, control_library_replicates,
                             ko_library_prefix, ko_library_replicates,
                             reference_genome_twobit, roi_file, cradle_out,   
                             control_flag=False):
    """Call peaks for control and ko library using cradle"""
    print("calling cradle peaks")
    input_rep_paths = list(map(lambda x: "_".join([input_library_prefix, x, "filtered"]) , input_library_replicates.split(" ")))
    control_rep_paths = list(map(lambda x: "_".join([control_library_prefix, x, "filtered"]) , control_library_replicates.split(" ")))
    ko_rep_paths = list(map(lambda x: "_".join([ko_library_prefix, x, "filtered"]) , ko_library_replicates.split(" ")))

    control_corrected_out = os.path.dirname(control_library_prefix)
    control_cradle_out = os.path.join(cradle_out, os.path.basename(control_library_prefix))

    ko_cradle_out = os.path.join(cradle_out, os.path.basename(ko_library_prefix))
    ko_corrected_out = os.path.dirname(ko_library_prefix)

    if control_flag:
        control_args = ["bash", f"{CURRENT_DIR_PATH}/call_peaks_cradle.sh"]
        control_args += [control_corrected_out, control_cradle_out]
        control_args += input_rep_paths
        control_args += control_rep_paths
        control_args += [roi_file, reference_genome_twobit]
        subprocess.call(control_args)

    ko_args = ["bash", f"{CURRENT_DIR_PATH}/call_peaks_cradle.sh"]
    ko_args = ["bash", f"{CURRENT_DIR_PATH}/call_peaks_cradle.sh"]
    ko_args += [ko_corrected_out, ko_cradle_out]
    ko_args += input_rep_paths
    ko_args += ko_rep_paths
    ko_args += [roi_file, reference_genome_twobit]
    subprocess.call(ko_args)
    return

# main script
def call_peaks(input_library_prefix, input_library_replicates, input_library_read_pairs, input_library_umi_index, input_library_suffix, 
               control_library_prefix, control_library_replicates, control_library_read_pairs, control_library_umi_index, control_library_suffix, 
               ko_library_prefix, ko_library_replicates, ko_library_read_pairs, ko_library_umi_index, ko_library_suffix,
               reference_genome, reference_genome_twobit, roi_file, 
               input_flag=False, control_flag=False, 
               dedup_flag=False, align_flag=False, bigwig_flag=False, filter_flag=False, 
               peak_flag=False, cradle_peak_flag=False):
    """Run everything from dedup to peak calling"""
    umi_flag = False
    if input_library_umi_index:
        umi_flag = True
        lib_suffix = "001.deduped.fastq"

        if dedup_flag:
            # use starrdust to remove duplicates
            remove_dups(input_library_prefix, input_library_replicates, input_library_read_pairs, input_library_umi_index, input_library_suffix, 
                        control_library_prefix, control_library_replicates, control_library_read_pairs, control_library_umi_index, control_library_suffix, 
                        ko_library_prefix, ko_library_replicates, ko_library_read_pairs, ko_library_umi_index, ko_library_suffix,
                        input_flag=input_flag, control_flag=control_flag)
        
        # overwrite library suffix
        input_library_suffix = lib_suffix
        control_library_suffix = lib_suffix
        ko_library_suffix = lib_suffix

    input_library_aligned_prefix = input_library_prefix.replace("raw_data", "aligned_reads")
    control_library_aligned_prefix = control_library_prefix.replace("raw_data", "aligned_reads")
    ko_library_aligned_prefix = ko_library_prefix.replace("raw_data", "aligned_reads")

    # align read files to the reference genome
    if align_flag:
        align_reads(input_library_prefix, input_library_replicates, input_library_read_pairs, input_library_suffix, input_library_aligned_prefix, 
                    control_library_prefix, control_library_replicates, control_library_read_pairs, control_library_suffix, control_library_aligned_prefix, 
                    ko_library_prefix, ko_library_replicates, ko_library_read_pairs, ko_library_suffix, ko_library_aligned_prefix, 
                    reference_genome, input_flag=input_flag, control_flag=control_flag)
    
    input_library_filtered_prefix = input_library_aligned_prefix.replace("aligned_reads", "filtered_libraries")
    control_library_filtered_prefix = control_library_aligned_prefix.replace("aligned_reads", "filtered_libraries")
    ko_library_filtered_prefix = ko_library_aligned_prefix.replace("aligned_reads", "filtered_libraries")

    # filter read files using samtools
    if filter_flag:
        filter_reads(input_library_aligned_prefix, input_library_replicates, input_library_filtered_prefix,
                     control_library_aligned_prefix, control_library_replicates, control_library_filtered_prefix,
                     ko_library_aligned_prefix, ko_library_replicates, ko_library_filtered_prefix, roi_file, 
                     umi=umi_flag, input_flag=input_flag, control_flag=control_flag)
    
    # create bigwig files from aligned and filtered bam files
    if bigwig_flag:
        create_bigwig(input_library_aligned_prefix, input_library_replicates,
                      control_library_aligned_prefix, control_library_replicates,
                      ko_library_aligned_prefix, ko_library_replicates,
                      input_flag=input_flag, control_flag=control_flag)

    # call peaks
    if peak_flag:
        call_peaks_helper(input_library_filtered_prefix, 
                          control_library_filtered_prefix, 
                          ko_library_filtered_prefix, 
                          control_flag=control_flag)

    if cradle_peak_flag:
        cradle_out_dir = os.path.dirname(input_library_aligned_prefix.replace("aligned_reads", "results/peaks"))
        os.makedirs(cradle_out_dir, exist_ok=True)
        call_peaks_helper_cradle(input_library_aligned_prefix, input_library_replicates,
                                 control_library_aligned_prefix, control_library_replicates,
                                 ko_library_aligned_prefix, ko_library_replicates,
                                 reference_genome_twobit, roi_file, cradle_out_dir,   
                                 control_flag=control_flag)        

    return 
