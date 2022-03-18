import os
import pandas as pd
import pybedtools
import pysam
import pyfastx

# TODO: take tmp dir from arguments
pybedtools.helpers.set_tempdir("/data5/deepro/tmp/")

####################
# Filename parsing #
####################

def get_prefix(filename):
    file_pre = os.path.splitext(filename)[0]
    return file_pre

# filtered bam filepaths
def get_filtered_bams(lib_prefix, lib_reps):
    """
    Function to fetch the filepaths of the filtered library replicates and merged bam file
    """
    lib_rep_pre = ["_".join([lib_prefix, r, "filtered"]) for r in lib_reps.split()]
    lib_rep_bams = list(map(lambda x: x.replace("raw_data", "aligned_reads")+".bam", lib_rep_pre))
    lib_merged_bam = lib_prefix.replace("raw_data", "filtered_libraries") + ".bam"
    return lib_rep_bams, lib_merged_bam

# coverage bed filepath
def get_coverage_files(rep_bams, merged_bam):
    """
    Function that outputs the filepath of where the coverage files will be stored based on the input bam filepath
    """
    rep_covs = [get_prefix(rb).replace("aligned_reads", "results/report/cov") + "_cov.bed" for rb in rep_bams]
    merged_cov = get_prefix(merged_bam).replace("filtered_libraries", "results/report/cov") + "_cov.bed"
    return rep_covs, merged_cov

###################
# Read QC Metrics #
###################

# roi info :: metric 1,2,3
def get_roi_info(roi_file):
    df_roi = pd.read_csv(roi_file, sep="\t")
    roi_num = len(df_roi)
    roi_total_bps = sum(df_roi.Size)
    roi_meansize = roi_total_bps//roi_num
    return roi_num, roi_meansize, roi_total_bps

# number of reads info :: metric 4,5
def get_num_reads_fastq(read_file):
    val = pysam.view("-c", read_file)    
    return val

def get_num_reads_bam(bam_file):
    val = pysam.view("-c", bamfile)
    return val

# overall coverage :: metric 8
def get_coverage(num_reads, len_region, read_length, pe=True):
    factor = 2 if pe else 1
    coverage = (num_reads*read_length*factor)/len_region
    return coverage

# roi individual coverage :: metric 9, 10
def get_roi_depth(filtered_bam, roi_sorted_bed, bed_out):
    bam = pybedtools.BedTool(filtered_bam)
    roi = pybedtools.BedTool(roi_sorted_bed)
    c = roi.coverage(bam)
    os.makedirs(os.path.dirname(bed_out), exist_ok=True)
    c.saveas(bed_out)
    return 
