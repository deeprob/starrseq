import os
import numpy as np
import pandas as pd
import pybedtools
import pysam
import multiprocessing as mp
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_theme(style="darkgrid")


# TODO: take tmp dir from arguments
pybedtools.helpers.set_tempdir("/data5/deepro/tmp/")

####################
# Filename parsing #
####################

def get_prefix(filename):
    file_pre = os.path.splitext(filename)[0]
    return file_pre

def get_html_file(out_pre, out_short):
    out_dir = os.path.dirname(out_pre).replace(f"raw_data/{out_short}", "results/report/htmls")
    os.makedirs(out_dir, exist_ok=True)
    out_html = os.path.join(out_dir, f"{out_short}.html")
    return out_html

# raw fastq filepaths
def get_raw_fastq(lib_prefix, lib_reps, lib_pairs, lib_suffix):
    """
    Function to fetch the filepaths of the paired end library replicates raw fastq files
    """
    lib_reps = ["_".join([lib_prefix, r, p, lib_suffix]) for r in lib_reps.split() for p in lib_pairs.split()]
    return [lib_reps[i] for i in range(0, len(lib_reps), 2)] if lib_pairs else lib_reps

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

def get_depth_figure_path(merged_cov_bed):
    return merged_cov_bed.replace(".bed", ".png")


###################
# Read QC Metrics #
###################

# roi info :: metric 1,2,3
def get_roi_info(roi_file):
    df_roi = pd.read_csv(roi_file, sep="\t", header=None)
    roi_num = len(df_roi)
    roi_total_bps = sum(df_roi.iloc[:,2].subtract(df_roi.iloc[:,1]))
    roi_meansize = roi_total_bps//roi_num
    return roi_num, roi_meansize, roi_total_bps

def add_roi_summary(roi_num, roi_meansize, roi_total_bps):
    roi_html = f"""
            <ol type="A" value="1">
            <li>The number of ROI is: <u>{roi_num}</u></li>
            <li>The mean size of ROI is: <u>{roi_meansize} bps</u></li> 
            <li>The total size of the ROI is: <u>{roi_total_bps} bps</u></li>
            </ol>
    """    
    return roi_html

# number of reads info :: metric 4,5
def get_num_reads_fastq(read_file):
    val = pysam.view("-c", read_file)    
    return int(val.strip())

def get_raw_reads_info(lib_prefix, lib_reps, lib_pairs, lib_suffix):
    lib_fastq_files = get_raw_fastq(lib_prefix, lib_reps, lib_pairs, lib_suffix)
    pool = mp.Pool(len(lib_fastq_files))
    lib_rep_reads = list(pool.map(get_num_reads_fastq, lib_fastq_files))
    pool.close()
    pool.join()
    return lib_rep_reads

def add_raw_reads_info(lib_rep_reads):
    raw_reads_html = "<ol type='A' value='1'>\n"
    for i,r in enumerate(lib_rep_reads):
        rep_num = i+1
        raw_reads_html += f"<li>Raw reads rep{rep_num}: <u>{r}</u></li>\n"
    raw_reads_html += "</ol>\n"
    return raw_reads_html

def get_num_reads_bam(bam_file):
    val = pysam.view("-c", bam_file)
    return int(val.strip())

def get_filtered_reads_info(lib_prefix, lib_reps):
    lib_bam_files, merged_bam = get_filtered_bams(lib_prefix, lib_reps)
    pool = mp.Pool(len(lib_bam_files) + 1)
    lib_bam_reads = list(pool.map(get_num_reads_bam, lib_bam_files + [merged_bam]))
    pool.close()
    pool.join()
    return lib_bam_reads

def add_filtered_reads_info(lib_bam_reads):
    filtered_reads_html = "<ol type='A' value='1'>\n"
    for i,r in enumerate(lib_bam_reads[:-1]):
        rep_num = i+1
        filtered_reads_html += f"<li>Filtered reads rep{rep_num}: <u>{r}</u></li>\n"
    filtered_reads_html += f"<li>Filtered reads merged: <u>{lib_bam_reads[-1]}</u></li>\n"
    filtered_reads_html += "</ol>\n"
    return filtered_reads_html

# overall coverage :: metric 8
def get_coverage(num_reads, len_region, read_length, pe=True):
    factor = 2 if pe else 1
    coverage = (num_reads*read_length*factor)/len_region
    return int(coverage)

def add_coverage_info(lib_coverage):
    coverage_html = "<ol type='A' value='1'>\n"
    for i,r in enumerate(lib_coverage[:-1]):
        rep_num = i+1
        coverage_html += f"<li>Coverage rep{rep_num}: <u>{r}x</u></li>\n"
    coverage_html += f"<li>Coverage merged: <u>{lib_coverage[-1]}x</u></li>\n"
    coverage_html += "</ol>\n"
    return coverage_html

# roi individual coverage :: metric 9, 10
def get_roi_depth(filtered_bam, roi_sorted_bed, bed_out):
    bam = pybedtools.BedTool(filtered_bam)
    roi = pybedtools.BedTool(roi_sorted_bed)
    c = roi.coverage(bam)
    os.makedirs(os.path.dirname(bed_out), exist_ok=True)
    c.saveas(bed_out)
    return

def make_coverage_bed(lib_prefix, lib_reps, roi):
    # lib bam files
    bam_reps, bam_merged = get_filtered_bams(lib_prefix, lib_reps)
    # input coverage beds
    rep_cov, merged_cov = get_coverage_files(bam_reps, bam_merged) 
    # make the coverage beds
    all_bams = bam_reps + [bam_merged]
    all_cov_beds = rep_cov + [merged_cov]
    map_iter = [(b,roi,c) for b,c in zip(all_bams, all_cov_beds)]
    pool = mp.Pool(8)
    pool.starmap(get_roi_depth, map_iter)
    pool.close()
    pool.join()
    return

def read_and_extract_coverage(cov_bed):
    df = pd.read_csv(cov_bed, sep="\t", header=None)
    roi_depth = df[3]
    return roi_depth

def atleast_num_rois(val, dfs):
    num_vals = (dfs>val).sum()
    return num_vals

def atleast_percent_rois(val, dfs):
    num_vals = (dfs>val).sum()
    percent_vals = round(num_vals*100/len(dfs), 3)
    return percent_vals

def plot_read_depth_helper(cov_beds, depth_fig):
    colnames = [f"rep{i}" for i in range(1, len(cov_beds))]
    colnames += ["merged"]

    dfs = pd.concat(list(map(read_and_extract_coverage, cov_beds)), axis=1)
    dfs.columns = colnames

    min_reads_vals = [i for i in range(10, 200, 10)]
    df = pd.DataFrame({"Minimum reads":min_reads_vals})

    df_roi_nums = pd.concat((df, df["Minimum reads"].apply(atleast_num_rois, args=(dfs, ))), axis=1).set_index("Minimum reads")
    df_roi_percent = pd.concat((df, df["Minimum reads"].apply(atleast_percent_rois, args=(dfs, ))), axis=1).set_index("Minimum reads")

    fig, axes = plt.subplots(1, 2, sharex=True, figsize=(15, 4))
    fig.suptitle("ROI Depth Statistics")
    sns.lineplot(data=df_roi_nums, linewidth=2, ax=axes[0])
    axes[0].set_ylabel("Number of ROIs")
    sns.lineplot(data=df_roi_percent, linewidth=2, ax=axes[1])
    axes[1].set_ylabel("Percentage of ROIs")
    plt.tight_layout()
    fig.savefig(depth_fig)
    return 

def plot_read_depth_stats(in_prefix, in_reps, out_prefix, out_reps, roi):
    """
    Plots the minimum number of reads aligned to ROIs vs Number of ROIs 
    """
    # make input coverage beds
    make_coverage_bed(in_prefix, in_reps, roi)
    # make output coverage beds
    make_coverage_bed(out_prefix, out_reps, roi)
    # input coverage bed and depth fig path
    input_bam_reps, input_bam_merged = get_filtered_bams(in_prefix, in_reps)
    input_rep_cov, input_merged_cov = get_coverage_files(input_bam_reps, input_bam_merged)
    in_cov = input_rep_cov + [input_merged_cov]
    in_depth_fig = get_depth_figure_path(input_merged_cov)
    # output coverage bed and depth fig path
    output_bam_reps, output_bam_merged = get_filtered_bams(out_prefix, out_reps)
    output_rep_cov, output_merged_cov = get_coverage_files(output_bam_reps, output_bam_merged)
    out_cov = output_rep_cov + [output_merged_cov]
    out_depth_fig = get_depth_figure_path(output_merged_cov)
    # input 
    plot_read_depth_helper(in_cov, in_depth_fig)
    # output
    plot_read_depth_helper(out_cov, out_depth_fig)
    return in_depth_fig , out_depth_fig

def add_depth_fig(depth_fig):
    depth_html = f"<img src={depth_fig} alt='cfg'>\n"
    return depth_html


######################
# Library QC Metrics #
######################

# within and between library correlation :: metrics 1,2,3

def plot_library_correlation_helper(all_cov_beds, fig_out):

    # get them as a dataframe
    df = pd.concat(list(map(read_and_extract_coverage, all_cov_beds)), axis=1)
    df.columns = [f"{line} Rep:{i}" for line in ["Input", "Output"] for i in range(1,len(all_cov_beds)//2+1)]

    ## create the subplot figure
    fig, axes = plt.subplots(
        nrows=len(all_cov_beds), ncols=len(all_cov_beds), 
        figsize=(12,8))

    for i in range(len(all_cov_beds)):
        for j in range(len(all_cov_beds)):
            if i==j:
                # diagonal
                lib_name = df.iloc[:,i].name
                axes[i][j].annotate(lib_name, xy=(0.5,0.5), xycoords='axes fraction', ha='center', fontsize=12, fontweight="bold")
                axes[i][j].set_axis_off()
            
            if i>j:
                # lower triangle
                sns.scatterplot(x=df.iloc[:,j], y= df.iloc[:,i], ax=axes[i][j])
                axes[i][j].set_ylabel(None)
                axes[i][j].set_xlabel(None)

                if j>0:
                    axes[i][j].yaxis.set_ticklabels([])


                if i<5:
                    axes[i][j].xaxis.set_ticklabels([])

            if i<j:
                # upper triangle
                r,p = pearsonr(df.iloc[:,i],df.iloc[:,j])
                axes[i][j].annotate('\u03C1 = {:.2f}'.format(r), xy=(0.5,0.5), xycoords='axes fraction', ha='center')
                axes[i][j].set_axis_off()


    plt.tight_layout()
    fig.savefig(fig_out)
    return fig_out

def get_corr_fig_path(out_pre, out_short):
    out_dir = os.path.dirname(out_pre).replace(f"raw_data/{out_short}", "results/report/corr")
    os.makedirs(out_dir, exist_ok=True)
    out_fig = os.path.join(out_dir, f"{out_short}.png")
    return out_fig

def plot_library_correlation(
    in_prefix, in_reps, 
    out_prefix, out_reps, out_short):
    # get input cov beds
    # lib bam files
    in_bam_reps, in_bam_merged = get_filtered_bams(in_prefix, in_reps)
    # input coverage beds
    in_rep_cov, in_merged_cov = get_coverage_files(in_bam_reps, in_bam_merged)
    # get output cov beds
    # lib bam files
    out_bam_reps, out_bam_merged = get_filtered_bams(out_prefix, out_reps)
    # input coverage beds
    out_rep_cov, out_merged_cov = get_coverage_files(out_bam_reps, out_bam_merged)
    all_rep_cov = in_rep_cov + out_rep_cov

    # get corr fig
    cfp = get_corr_fig_path(out_prefix, out_short)
    cfp = plot_library_correlation_helper(all_rep_cov, cfp)
    return cfp

def add_correlation_fig(corr_fig):
    corr_html = f"<img src={corr_fig} alt='cfg'>\n"
    return corr_html

# fold change correlation between replicates :: metrics 4

def get_repnum(in_beds, ko_beds):
    assert len(in_beds) == len(ko_beds)
    return len(in_beds)

def read_coverage_convert_to_rpkm(cov_bed):
    """
    Reads the location and read depth of a region from a bed file
    Normalizes the read depth to RPKM value
    """
    df = pd.read_csv(cov_bed, sep="\t", header=None, usecols=[0,1,2,3], names=["chrm", "start", "end", "reads"])
    # assign pseudo count of 1 to 0 read regions
    df.reads = df.reads.replace(0, 1)
    # get the length of the regions to calculate rpkm values
    df_gene_length = df.end - df.start
    # calculate rpkm
    df_norm_reads = (df.reads*10**6*10**3)/(df_gene_length*df.reads.sum())
    df["rpkm"] = df_norm_reads
    df = df.set_index(["chrm", "start", "end"])
    return df.loc[:, ["rpkm"]]

def convert_beds_to_rpkmdf(in_beds, out_beds):
    rep_num = get_repnum(in_beds, out_beds)
    df = pd.concat(list(map(read_coverage_convert_to_rpkm, in_beds+out_beds)), axis=1)
    df.columns = [f"{line} Rep:{i}" for line in ["Input", "Output"] for i in range(1, rep_num + 1)]
    return df

def get_log2_rpkm_fc(df, in_col, out_col):
    df_rep = df.iloc[:, [in_col, out_col]]
    df_lfc = np.log2(df_rep.iloc[:,1]/df_rep.iloc[:,0])
    return df_lfc

def get_lfc_df(df):
    rep_num = df.shape[1]//2
    df_lfc = pd.concat(
        [get_log2_rpkm_fc(df, i, i + rep_num) for i in range(rep_num)], 
        axis=1
        )
    df_lfc.columns = [f"Rep {i}" for i in range(1, rep_num + 1)]
    return df_lfc

def get_fc_fig_path(out_pre, out_short):
    out_dir = os.path.dirname(out_pre).replace(f"raw_data/{out_short}", "results/report/fc")
    os.makedirs(out_dir, exist_ok=True)
    out_fig = os.path.join(out_dir, f"{out_short}.png")
    return out_fig

def plot_fc_correlation_helper(in_beds, out_beds, fig_out):

    # get them as a dataframe
    df = convert_beds_to_rpkmdf(in_beds, out_beds)
    lfc_df = get_lfc_df(df)

    ## create the subplot figure
    fig, axes = plt.subplots(
        nrows=len(in_beds), ncols=len(in_beds), 
        figsize=(8,6))

    for i in range(len(in_beds)):
        for j in range(len(in_beds)):
            if i==j:
                # diagonal
                lib_name = lfc_df.iloc[:,i].name
                axes[i][j].annotate(lib_name, xy=(0.5,0.5), xycoords='axes fraction', ha='center', fontsize=12, fontweight="bold")
                axes[i][j].set_axis_off()
            
            if i>j:
                # lower triangle
                sns.scatterplot(x=lfc_df.iloc[:,j], y= lfc_df.iloc[:,i], ax=axes[i][j])
                axes[i][j].set_ylabel(None)
                axes[i][j].set_xlabel(None)

                if j>0:
                    axes[i][j].yaxis.set_ticklabels([])


                if i<len(in_beds)-1:
                    axes[i][j].xaxis.set_ticklabels([])

            if i<j:
                # upper triangle
                r,p = pearsonr(lfc_df.iloc[:,i],lfc_df.iloc[:,j])
                axes[i][j].annotate('\u03C1 = {:.2f}'.format(r), xy=(0.5,0.5), xycoords='axes fraction', ha='center')
                axes[i][j].set_axis_off()


    plt.tight_layout()
    fig.savefig(fig_out)
    return fig_out

def plot_fc_correlation(
    in_prefix, in_reps, 
    out_prefix, out_reps, out_short):
    # get input cov beds
    # lib bam files
    in_bam_reps, in_bam_merged = get_filtered_bams(in_prefix, in_reps)
    # input coverage beds
    in_rep_cov, in_merged_cov = get_coverage_files(in_bam_reps, in_bam_merged)
    # get output cov beds
    # lib bam files
    out_bam_reps, out_bam_merged = get_filtered_bams(out_prefix, out_reps)
    # input coverage beds
    out_rep_cov, out_merged_cov = get_coverage_files(out_bam_reps, out_bam_merged)

    # get corr fig
    cfp = get_fc_fig_path(out_prefix, out_short)
    cfp = plot_fc_correlation_helper(in_rep_cov, out_rep_cov, cfp)
    return cfp

def add_fc_fig(fc_fig):
    fc_html = f"<img src={fc_fig} alt='cfg'>\n"
    return fc_html
