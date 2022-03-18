import os
import multiprocessing as mp
import quality_control.utils as ut


def make_coverage_bed(in_prefix, in_reps, out_prefix, out_reps, roi):
    # input bam files
    input_bam_reps, input_bam_merged = ut.get_filtered_bams(in_prefix, in_reps)
    # output bam files
    output_bam_reps, output_bam_merged = ut.get_filtered_bams(out_prefix, out_reps)
    # input coverage beds
    input_rep_cov, input_merged_cov = ut.get_coverage_files(input_bam_reps, input_bam_merged) 
    # output coverage beds
    output_rep_cov, output_merged_cov = ut.get_coverage_files(output_bam_reps, output_bam_merged)
    all_bams = input_bam_reps + [input_bam_merged] + output_bam_reps + [output_bam_merged]
    all_cov_beds = input_rep_cov + [input_merged_cov] + output_rep_cov + [output_merged_cov]
    map_iter = [(b,roi,c) for b,c in zip(all_bams, all_cov_beds)]
    pool = mp.Pool(8)
    pool.starmap(ut.get_roi_depth, map_iter)
    pool.close()
    pool.join()
    return


def generate_read_qc_report(in_prefix, in_reps, out_prefix, out_reps, roi):
    print("Generating Report ...")
    # TODO: first write all roi information

    # TODO: then write all reads information

    print("Making ROI coverage beds ...")
    make_coverage_bed(in_prefix, in_reps, out_prefix, out_reps, roi)

    # TODO: now write all roi coverage information

    return


def generate_report(in_prefix, in_reps, out_prefix, out_reps, roi):

    generate_read_qc_report(in_prefix, in_reps, out_prefix, out_reps, roi)

    return
