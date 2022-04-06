# Import modules
import argparse
import os

# internal imports
import reads_to_peaks.utils as ut


CURRENT_DIR_PATH = os.path.dirname(os.path.abspath(__file__))


def remove_dups(input_library_prefix, input_library_replicates, input_library_read_pairs, input_library_umi_index, input_library_suffix,
                output_library_prefix, output_library_replicates, output_library_read_pairs, output_library_umi_index, output_library_suffix,
                input_flag, output_flag):
    """
    Remove duplicates from all three libraries using starrdust
    """
    if input_flag:
        ut.get_deduped_files(input_library_prefix, input_library_replicates, input_library_read_pairs, input_library_umi_index, input_library_suffix)
    if output_flag:
        ut.get_deduped_files(output_library_prefix, output_library_replicates, output_library_read_pairs, output_library_umi_index, output_library_suffix)
    return


def align_reads(input_library_prefix, input_library_replicates, input_library_read_pairs, input_library_suffix, input_library_aligned_prefix,
                output_library_prefix, output_library_replicates, output_library_read_pairs, output_library_suffix, output_library_aligned_prefix,
                reference_genome, input_flag, output_flag, cores=64):
    """Align reads to the reference genome"""

    if input_flag:
        # align input reads to the reference genome
        ut.align_reads_helper(reference_genome, input_library_prefix, input_library_replicates, input_library_read_pairs, input_library_suffix, input_library_aligned_prefix, cores)
    if output_flag:
        # align output reads to the reference genome
        ut.align_reads_helper(reference_genome, output_library_prefix, output_library_replicates, output_library_read_pairs, output_library_suffix, output_library_aligned_prefix, cores)
    return


def filter_reads(input_library_aligned_prefix, input_library_replicates, input_library_filtered_prefix,
                output_library_aligned_prefix, output_library_replicates, output_library_filtered_prefix,
                roi_file, 
                umi, input_flag, output_flag):
    """filter bad reads"""
    dedup = "true"
    if umi:
        dedup = "false"

    if input_flag:
        # filter input reads samtools -F 2828 or 2852 (depending on the use of starrdust or picard), -f 2 -q 30
        ut.filter_reads_helper(input_library_aligned_prefix, input_library_replicates, input_library_filtered_prefix, roi_file, dedup)
    if output_flag:
        # filter output reads
        ut.filter_reads_helper(output_library_aligned_prefix, output_library_replicates, output_library_filtered_prefix, roi_file, dedup)
    return 


def call_peaks(input_library_filtered_prefix, output_library_filtered_prefix, 
            input_library_aligned_prefix, input_library_replicates,
            output_library_aligned_prefix, output_library_replicates,
            reference_genome, reference_genome_twobit, roi_file,
            starrpeaker_data_dir, cradle_data_dir,
            input_flag, output_flag, 
            starrpeaker_peak_flag, cradle_peak_flag, macs2_peak_flag, replicate_peak_flag):

    if starrpeaker_peak_flag:
        ut.call_starrpeaker_peaks(input_library_filtered_prefix, 
            output_library_filtered_prefix,
            starrpeaker_data_dir,
            output_flag=output_flag)

    if cradle_peak_flag:
        # create bigwig files from aligned and filtered bam files, required for cradle
        ut.create_bigwig(input_library_aligned_prefix, input_library_replicates,
            output_library_aligned_prefix, output_library_replicates,
            input_flag=input_flag, output_flag=output_flag)

        ut.call_cradle_peaks(input_library_aligned_prefix, input_library_replicates,
                            output_library_aligned_prefix, output_library_replicates,
                            reference_genome_twobit, roi_file, cradle_data_dir,
                            output_flag=output_flag)

    if macs2_peak_flag:
        ut.call_macs2_peaks(input_library_filtered_prefix, 
            output_library_filtered_prefix,
            output_flag=output_flag)
    
    #TODO: peak call per replicate for starrpeaker peak caller to compare peaks called between replicates
    if replicate_peak_flag:
        ut.call_starrpeaker_peaks_for_each_replicate(
            input_library_aligned_prefix, input_library_replicates, 
            output_library_aligned_prefix, output_library_replicates,
            starrpeaker_data_dir,
            output_flag=output_flag)        

    return 

def reads2peaks(input_library_prefix, input_library_replicates, input_library_read_pairs, input_library_umi_index, input_library_suffix, 
               output_library_prefix, output_library_replicates, output_library_read_pairs, output_library_umi_index, output_library_suffix,
               reference_genome, reference_genome_twobit, roi_file,
               starrpeaker_data_dir, cradle_data_dir, 
               input_flag=True, output_flag=True, 
               dedup_flag=True, align_flag=True, filter_flag=True, 
               starrpeaker_peak_flag=True, cradle_peak_flag=True, macs2_peak_flag=True,
               rep_peak_flag=True):
    """
    Run everything from dedup to peak calling
    """

    # deduplication
    umi_flag = False
    if input_library_umi_index:
        umi_flag = True
        lib_suffix = "001.deduped.fastq"

        if dedup_flag:
            # use starrdust to remove duplicates
            remove_dups(input_library_prefix, input_library_replicates, input_library_read_pairs, input_library_umi_index, input_library_suffix, 
                        output_library_prefix, output_library_replicates, output_library_read_pairs, output_library_umi_index, output_library_suffix,
                        input_flag=input_flag, output_flag=output_flag)
        
        # overwrite library suffix
        input_library_suffix = lib_suffix
        output_library_suffix = lib_suffix

    input_library_aligned_prefix = input_library_prefix.replace("raw_data", "aligned_reads")
    output_library_aligned_prefix = output_library_prefix.replace("raw_data", "aligned_reads")

    # align
    if align_flag:
        align_reads(input_library_prefix, input_library_replicates, input_library_read_pairs, input_library_suffix, input_library_aligned_prefix, 
                    output_library_prefix, output_library_replicates, output_library_read_pairs, output_library_suffix, output_library_aligned_prefix,
                    reference_genome, input_flag=input_flag, output_flag=output_flag)
    
    input_library_filtered_prefix = input_library_aligned_prefix.replace("aligned_reads", "filtered_libraries")
    output_library_filtered_prefix = output_library_aligned_prefix.replace("aligned_reads", "filtered_libraries")

    # filter
    if filter_flag:
        filter_reads(input_library_aligned_prefix, input_library_replicates, input_library_filtered_prefix,
                     output_library_aligned_prefix, output_library_replicates, output_library_filtered_prefix,
                     roi_file, 
                     umi=umi_flag, input_flag=input_flag, output_flag=output_flag)

    # call peaks
    call_peaks(input_library_filtered_prefix, output_library_filtered_prefix, 
                input_library_aligned_prefix, input_library_replicates,
                output_library_aligned_prefix, output_library_replicates,
                reference_genome, reference_genome_twobit, roi_file,
                starrpeaker_data_dir, cradle_data_dir, 
                input_flag, output_flag, 
                starrpeaker_peak_flag, cradle_peak_flag, macs2_peak_flag, rep_peak_flag)

    return 

