import argparse
import prepare_meta.create_metafile as cm
import reads_to_peaks.utils as rpu
import reads_to_peaks.reads_to_peaks as rpp
import quality_control.report as qc


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="STARRSeq analysis")
    # only required argument, the meta file that stores library information
    parser.add_argument("meta_file", type=str, help="The meta json filepath where library information is stored, look at $prepare meta$ package")
    # calling the k or ko argument will start the reads to peaks process of the pipeline
    parser.add_argument("-k", "--ko", type=str, help="The output library name as given in the meta file, if present, this will run reads to peaks")
    parser.add_argument("-i", "--input", action="store_false", help="Do not run the input library pipeline")
    parser.add_argument("-o", "--output", action="store_false", help="Do not run the output library pipeline")
    parser.add_argument("-d", "--dedup", action="store_false", help="Do not run the deduplication pipeline")
    parser.add_argument("-a", "--align", action="store_false", help="Do not run the alignment pipeline")
    parser.add_argument("-f", "--filter", action="store_false", help="Do not run the filter pipeline")
    parser.add_argument("-s", "--starrpeaker", action="store_false", help="Do not run the starrpeaker peak calling pipeline")
    parser.add_argument("-c", "--cradle", action="store_false", help="Do not run the cradle peak calling pipeline")
    parser.add_argument("-m", "--macs2", action="store_false", help="Do not run the macs2 peak calling pipeline")
    # calling q or qual argument will start the quality control process of the pipeline
    parser.add_argument("-q", "--qual", action="store_true", help="Generate a quality control report of the specified input and ko library")

    cli_args = parser.parse_args()

    # deciding which sub-package to run
    if cli_args.ko:
        # generating the arguments
        args = rpu.create_args(cli_args.meta_file, cli_args.ko, 
                    cli_args.input, cli_args.output, 
                    cli_args.dedup, cli_args.align, cli_args.filter, 
                    cli_args.starrpeaker, cli_args.cradle, cli_args.macs2)

        if cli_args.qual:
            # running the quality control pipeline with the qual flag specified
            qc.generate_report(
                args.input_library_prefix, args.input_library_reps, 
                args.output_library_prefix, args.output_library_reps,
                args.roi_file)

        else:
            # running the reads to peaks pipeline when only the ko line is specified
            rpp.reads2peaks(
                args.input_library_prefix, args.input_library_reps, args.input_library_pair, args.input_library_umi, args.input_library_suffix, 
                args.output_library_prefix, args.output_library_reps, args.output_library_pair, args.output_library_umi, args.output_library_suffix,
                args.reference_genome, args.reference_genome_twobit, args.roi_file,
                args.starrpeaker_data_dir, args.cradle_data_dir, 
                args.input_flag, args.output_flag, 
                args.dedup_flag, args.align_flag, args.filter_flag, 
                args.starrpeaker_peak_flag, args.cradle_peak_flag, args.macs2_peak_flag)

    else:
        # with no flag specified, create a meta file that fills a dictionary with path to starrseq libraries
        cm.create_meta(cli_args.meta_file, force=True)
