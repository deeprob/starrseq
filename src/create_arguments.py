import os
import json
from argparse import Namespace


def create_args(meta_file, ko_name, ref_gen, roi_master, roi_sorted, run=False, input_flag=True, control_flag=True):

    with open(meta_file, "r") as f: 
        meta_dict = json.load(f)

    args = Namespace(
        # from metadata and globals
        input_library_prefix = meta_dict["input"]["prefix"],
        control_library_prefix = meta_dict["control"]["prefix"],
        ko_library_prefix = meta_dict[ko_name]["prefix"],
        input_library_reps = meta_dict["input"]["replicates"],
        control_library_reps = meta_dict["control"]["replicates"],
        ko_library_reps = meta_dict[ko_name]["replicates"],
        input_library_pair = meta_dict["input"]["read_pairs"],
        control_library_pair= meta_dict["control"]["read_pairs"],
        ko_library_pair = meta_dict[ko_name]["read_pairs"],
        input_library_umi = meta_dict["input"]["umi"],
        control_library_umi = meta_dict["control"]["umi"],
        ko_library_umi = meta_dict[ko_name]["umi"],
        input_library_suffix = meta_dict["input"]["suffix"],
        control_library_suffix = meta_dict["control"]["suffix"],
        ko_library_suffix = meta_dict[ko_name]["suffix"],
        reference_genome = ref_gen,
        region_of_interest_master = roi_master,
        region_of_interest_sorted = roi_sorted,
        ko_name = "eGFP-ATF2", # TODO: add this information in metadata and get this information from metadata

        # for internal use
        region_of_interest_sorted_seq = roi_sorted.replace("sorted", "sorted.seq"),
        region_of_interest_homer_background = roi_sorted.replace("sorted", "homer"),

        ## align and filter
        input_library_aligned_prefix = meta_dict["input"]["prefix"].replace("raw_data", "aligned_reads"),
        control_library_aligned_prefix = meta_dict["control"]["prefix"].replace("raw_data", "aligned_reads"),
        ko_library_aligned_prefix = meta_dict[ko_name]["prefix"].replace("raw_data", "aligned_reads"),
        input_library_filtered_prefix = meta_dict["input"]["prefix"].replace("raw_data", "filtered_libraries"),
        control_library_filtered_prefix = meta_dict["control"]["prefix"].replace("raw_data", "filtered_libraries"),
        ko_library_filtered_prefix = meta_dict[ko_name]["prefix"].replace("raw_data", "filtered_libraries"),
        input_library_filtered_bam = meta_dict["input"]["prefix"].replace("raw_data", "filtered_libraries") + ".bam",
        control_library_filtered_bam = meta_dict["control"]["prefix"].replace("raw_data", "filtered_libraries") + ".bam",
        ko_library_filtered_bam = meta_dict[ko_name]["prefix"].replace("raw_data", "filtered_libraries") + ".bam",
        input_library_coverage_bed = meta_dict["input"]["prefix"].replace("raw_data", "filtered_libraries") + ".coverage.bed",
        control_library_coverage_bed = meta_dict["control"]["prefix"].replace("raw_data", "filtered_libraries") + ".coverage.bed",
        ko_library_coverage_bed = meta_dict[ko_name]["prefix"].replace("raw_data", "filtered_libraries") + ".coverage.bed",
        input_library_depth_bed = meta_dict["input"]["prefix"].replace("raw_data", "filtered_libraries") + ".depth.bed",
        control_library_depth_bed = meta_dict["control"]["prefix"].replace("raw_data", "filtered_libraries") + ".depth.bed",
        ko_library_depth_bed = meta_dict[ko_name]["prefix"].replace("raw_data", "filtered_libraries") + ".depth.bed",
        ## starrpeaker peaks
        control_peaks_prefix = os.path.join(meta_dict["control"]["prefix"].replace("raw_data", "results/peaks"), "starrpeaker"), 
        ko_peaks_prefix = os.path.join(meta_dict[ko_name]["prefix"].replace("raw_data", "results/peaks"), "starrpeaker"),
        control_peak_file = os.path.join(meta_dict["control"]["prefix"].replace("raw_data", "results/peaks"), "starrpeaker") + ".peak.final.bed",
        ko_peak_file = os.path.join(meta_dict[ko_name]["prefix"].replace("raw_data", "results/peaks"), "starrpeaker") + ".peak.final.bed",
        ko_dactive_file = os.path.join(meta_dict[ko_name]["prefix"].replace("raw_data", "results/peaks"), "starrpeaker") + ".peak.direct.active.bed",
        ko_dinactive_file = os.path.join(meta_dict[ko_name]["prefix"].replace("raw_data", "results/peaks"), "starrpeaker") + ".peak.direct.inactivate.bed",
        ko_daannotated_file = os.path.join(meta_dict[ko_name]["prefix"].replace("raw_data", "results/peaks"), "starrpeaker") + ".peak.dactive.annotated.bed",
        ko_diannotated_file = os.path.join(meta_dict[ko_name]["prefix"].replace("raw_data", "results/peaks"), "starrpeaker") + ".peak.dinactive.annotated.bed",
        ## cradle peaks
        cradle_control_peaks_prefix = meta_dict["control"]["prefix"].replace("raw_data", "results/peaks"),
        cradle_ko_peaks_prefix = meta_dict[ko_name]["prefix"].replace("raw_data", "results/peaks"),
        cradle_control_peak_file = os.path.join(meta_dict["control"]["prefix"].replace("raw_data", "results/peaks"), "CRADLE_peaks"),
        cradle_ko_peak_file = os.path.join(meta_dict[ko_name]["prefix"].replace("raw_data", "results/peaks"), "CRADLE_peaks"),
        cradle_control_activated_file = os.path.join(meta_dict["control"]["prefix"].replace("raw_data", "results/peaks"), "CRADLE_peaks") + "_activated",
        cradle_control_repressed_file = os.path.join(meta_dict["control"]["prefix"].replace("raw_data", "results/peaks"), "CRADLE_peaks") + "_repressed",
        cradle_ko_activated_file = os.path.join(meta_dict[ko_name]["prefix"].replace("raw_data", "results/peaks"), "CRADLE_peaks") + "_activated",
        cradle_ko_repressed_file = os.path.join(meta_dict[ko_name]["prefix"].replace("raw_data", "results/peaks"), "CRADLE_peaks") + "_repressed",
        cradle_ko_dactive_file = os.path.join(meta_dict[ko_name]["prefix"].replace("raw_data", "results/peaks"), "CRADLE_peaks") + "_direct_active",
        cradle_ko_dinactive_file = os.path.join(meta_dict[ko_name]["prefix"].replace("raw_data", "results/peaks"), "CRADLE_peaks") + "_direct_inactivate",
        cradle_ko_daannotated_file = os.path.join(meta_dict[ko_name]["prefix"].replace("raw_data", "results/peaks"), "CRADLE_peaks") + "_dactive_annotated",
        cradle_ko_diannotated_file = os.path.join(meta_dict[ko_name]["prefix"].replace("raw_data", "results/peaks"), "CRADLE_peaks") + "_dinactive_annotated",
        ## starrpeaker mea
        ko_dactive_mea_dir = os.path.join(meta_dict[ko_name]["prefix"].replace("raw_data", "results/mea"), "starrpeaker/active"),
        ko_dinactive_mea_dir = os.path.join(meta_dict[ko_name]["prefix"].replace("raw_data", "results/mea"), "starrpeaker/inactive"),
        ##cradle mea
        cradle_ko_dactive_mea_dir = os.path.join(meta_dict[ko_name]["prefix"].replace("raw_data", "results/mea"), "cradle/active"),
        cradle_ko_dinactive_mea_dir = os.path.join(meta_dict[ko_name]["prefix"].replace("raw_data", "results/mea"), "cradle/inactive"),

        # from global flags
        umi_flag = True if meta_dict["input"]["umi"] else False,
        input_flag = input_flag,
        control_flag = control_flag,
        run = run,

    )
    return args

