import os
import multiprocessing as mp
import quality_control.utils as ut
from itertools import starmap


def read_qc_html(
    roi_num, roi_meansize, roi_total_bps, 
    in_raw_reads, in_bam_reads, in_coverage, in_depth_fig,
    out_raw_reads, out_bam_reads, out_coverage, out_depth_fig,
    ):

    roi_html = ut.add_roi_summary(roi_num, roi_meansize, roi_total_bps)
    in_rr_html = ut.add_raw_reads_info(in_raw_reads)
    out_rr_html = ut.add_raw_reads_info(out_raw_reads)
    in_fr_html = ut.add_filtered_reads_info(in_bam_reads)
    out_fr_html = ut.add_filtered_reads_info(out_bam_reads)
    in_cv_html = ut.add_coverage_info(in_coverage)
    out_cv_html = ut.add_coverage_info(out_coverage)
    in_dp_html = ut.add_depth_fig(in_depth_fig)
    out_dp_html = ut.add_depth_fig(out_depth_fig)

    html_out = f"""
        <h1>Read Quality Control</h1>
        <h2>Regions of Interest (ROI) Summary</h2>
        {roi_html}
        <h2>Raw Reads Statistics</h2>
        <h3>Input</h3>
        {in_rr_html}
        <h3>Output</h3>
        {out_rr_html}
        <h2>Filtered Reads Statistics</h2>
        <h3>Input</h3>
        {in_fr_html}
        <h3>Output</h3>
        {out_fr_html}
        <h2>Read Depth Statistics</h2>
        <h3>Input coverage</h3>
        {in_cv_html}
        <h3>Minimum Input Reads per ROI</h3>
        {in_dp_html}
        <h3>Output coverage</h3>
        {out_cv_html}
        <h3>Minimum Output Reads per ROI</h3>
        {out_dp_html}          
        """
    return html_out



def generate_read_qc_report(
    in_prefix, in_reps, in_pairs, in_suffix, 
    out_prefix, out_reps, out_pairs, out_suffix, 
    roi):
    print("Generating Read QC Report ...")
    # ROI information
    roi_num, roi_meansize, roi_total_bps = ut.get_roi_info(roi)
    # raw reads info
    in_raw_reads = ut.get_raw_reads_info(in_prefix, in_reps, in_pairs, in_suffix)
    out_raw_reads = ut.get_raw_reads_info(out_prefix, out_reps, out_pairs, out_suffix)
    # filtered reads info
    in_bam_reads = ut.get_filtered_reads_info(in_prefix, in_reps)
    out_bam_reads = ut.get_filtered_reads_info(out_prefix, out_reps)
    # coverage statistics
    in_iter_map = [(ibr, roi_total_bps, 150) for ibr in in_bam_reads]
    in_coverage = list(starmap(ut.get_coverage, in_iter_map))
    out_iter_map = [(obr, roi_total_bps, 150) for obr in out_bam_reads]
    out_coverage = list(starmap(ut.get_coverage, out_iter_map))
    # depth statistics
    in_depth_fig, out_depth_fig = ut.plot_read_depth_stats(
        in_prefix, in_reps, 
        out_prefix, out_reps, 
        roi)

    rqc_html = read_qc_html( 
        roi_num, roi_meansize, roi_total_bps, 
        in_raw_reads, in_bam_reads, in_coverage, in_depth_fig,
        out_raw_reads, out_bam_reads, out_coverage, out_depth_fig,
        )


    return rqc_html


def lib_qc_html(
    corr_fig_path, fc_fig_path, 
    ):

    corr_html = ut.add_correlation_fig(corr_fig_path)
    fc_html = ut.add_fc_fig(fc_fig_path)

    html_out = f"""
        <h1>Library Quality Control</h1>
        <h2>Replicate Correlation Plots</h2>
        {corr_html}
        <h2>Fold Change Correlation</h2>
        {fc_html}
        <h2>Enhancer Summit Correlation</h2>

        <h2>Reproducible Peaks Statistics</h2>   
        """
    return html_out

def generate_lib_qc_report(
    in_prefix, in_reps, 
    out_prefix, out_reps, out_short):

    print("Generating Library QC Report ...")
    corr_fig_path = ut.plot_library_correlation(
        in_prefix, in_reps, 
        out_prefix, out_reps, out_short)

    fc_fig_path = ut.plot_fc_correlation(
        in_prefix, in_reps, 
        out_prefix, out_reps, out_short)

    lqc_html = lib_qc_html(corr_fig_path, fc_fig_path)


    return lqc_html


def generate_report(
    out_short, 
    in_prefix, in_reps, in_pairs, in_suffix, 
    out_prefix, out_reps, out_pairs, out_suffix, 
    roi):

    html_file = ut.get_html_file(out_prefix, out_short)

    read_qc_html = generate_read_qc_report(
        in_prefix, in_reps, in_pairs, in_suffix, 
        out_prefix, out_reps, out_pairs, out_suffix, roi)

    lib_qc_html = generate_lib_qc_report(
        in_prefix, in_reps,
        out_prefix, out_reps, out_short)


    # create html file path using short form info from meta dict
    with open(html_file, "w") as htmlf:
        htmlf.write(f"""
        <html>
        <header style="font-size: 30px;">
        STARRSeq Quality Control Report
        </header>
        
        <body>

        {read_qc_html}

        {lib_qc_html}

        </body>

        </html>
        """)

    return
