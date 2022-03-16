import os
import re
import subprocess
import argparse
import pybedtools
import pandas as pd
import seaborn as sns
import multiprocessing as mp

# take tmp dir from arguments
pybedtools.helpers.set_tempdir("/data5/deepro/tmp/")

CURRENT_DIR_PATH = os.path.dirname(os.path.abspath(__file__))

def read_tsv(filename):
    filebase = os.path.basename(filename)
    pattern = "(.+)\.window\.coverage\.bed"
    m = re.match(pattern, filebase)
    counts_colname = m.group(1)
    df = pd.read_csv(filename, sep="\t", header=None, usecols=[0,1,3])
    df.columns = ["chr", "start", counts_colname]
    df["end"] = df["start"] + 500
    df.set_index(["chr", "start", "end"], inplace=True)
    return df

def make_enhancer_fragments(roi_file, window_file, window_size=500, window_stride=50):
    """
    Break the ROIs into fragments of an user defined window size and stride
    """
    bin = pybedtools.BedTool().window_maker(b=roi_file ,w=window_size, s=window_stride)
    bin.saveas(window_file)
    return


def get_roi_coverage(file_a, file_b, file_out):
    """
    Function to get the coverage of the Enhancer Fragments
    file_a: ROI window file
    file_b: filtered bam file
    file_c: coverage file
    """
    exit_code = subprocess.call(["bash", f"{CURRENT_DIR_PATH}/get_coverage.sh", f"{file_a}", f"{file_b}", f"{file_out}", "false"])
    assert exit_code==0
    return


def main(roi_file, bam_files, wsize, wstride, force_flag):

    # Step 1: Make the window file
    window_file = os.path.splitext(roi_file)[0] + ".window.bed"
    if not os.path.exists(window_file) or force_flag:
        make_enhancer_fragments(roi_file, window_file, wsize, wstride)
    else:
        print(f"Window file {window_file} already present")
    # Step 2: Get the coverage of the window file against the bam file(s)
    coverage_files = [os.path.splitext(bam)[0] + ".window.coverage.bed" for bam in bam_files]
    cores = len(bam_files)
    pool = mp.Pool(cores)
    args_map = [(window_file, bam, cov) for bam,cov in zip(bam_files, coverage_files)]
    pool.starmap(get_roi_coverage, args_map)
    pool.close()
    pool.join()
    # Step 3: Store the coverage info for all bamfiles as a csv file
    mega_counts_file = os.path.join(os.path.dirname(bam_files[0]), "mega_counts.csv")
    dfs = list(map(read_tsv, coverage_files))
    mega_df = pd.concat(dfs, axis=1)
    mega_df.to_csv(mega_counts_file)
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create heatmap")

    parser.add_argument("roi", type=str, help="The ROI enhancer file")
    parser.add_argument("bam", type=str, help="Bamfile(s) to get the reads of the ROI windows", nargs="+")
    parser.add_argument("--size", type=int, help="Window size for the fragments", default=500)
    parser.add_argument("--stride", type=int, help="Window stride for the fragments", default=50)
    # TODO: mega counts file output
    parser.add_argument("-f", "--force", action="store_true")

    args = parser.parse_args()

    main(args.roi, args.bam, args.size, args.stride, args.force) 

    # example command
    # python window_coverage.py /data5/deepro/starrseq/computational_pipeline/data/master.sorted.bed /data5/deepro/starrseq/main_lib/filtered_libraries/16P12_1.bam /data5/deepro/starrseq/main_lib/filtered_libraries/ATF2.bam /data5/deepro/starrseq/main_lib/filtered_libraries/CC.bam /data5/deepro/starrseq/main_lib/filtered_libraries/CTCF.bam /data5/deepro/starrseq/main_lib/filtered_libraries/FOXA1.bam /data5/deepro/starrseq/main_lib/filtered_libraries/Input_SeqReady.bam /data5/deepro/starrseq/main_lib/filtered_libraries/LEF1.bam /data5/deepro/starrseq/main_lib/filtered_libraries/SCRT1.bam /data5/deepro/starrseq/main_lib/filtered_libraries/TCF7L2.bam 
