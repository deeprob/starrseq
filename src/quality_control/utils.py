import pandas as pd
import pybedtools
import pysam
import pyfastx


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
def 