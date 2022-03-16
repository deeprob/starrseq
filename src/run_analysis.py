import argparse
import reads_to_peaks as rp

def main():

    # check if umi indices exist
    if library_umis:
        input_umi = input_library_umi_index
        control_umi = control_library_umi_index
        ko_umi = ko_library_umi_index
    else:
        input_umi = ""
        control_umi = ""
        ko_umi = ""

    # call peaks
    rp.call_peaks(input_lp, input_lr, input_lrp, input_umi, input_ls, 
                  control_lp, control_lr, control_lrp, control_umi, control_ls,  
                  ko_lp, ko_lr, ko_lrp, ko_umi, ko_ls, 
                  reference_genome, reference_genome_twobit, roi_file,
                  input_flag=finput, control_flag=fcontrol, 
                  dedup_flag=fdedup, align_flag=falign, bigwig_flag=fbigwig, filter_flag=ffilter, 
                  peak_flag=fpeak, cradle_peak_flag=fcpeak)
    
    return 


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=
    "STARRSeq analysis;\
    Read files are expected to have the following structure:\
    {lp}_{replicate}_{read_pair}_{ls}")

    parser.add_argument("--input_prefix", type=str, help="prefix for the input library read file path")
    parser.add_argument("--control_prefix", type=str, help="prefix for the control library read file path")
    parser.add_argument("--ko_prefix", type=str, help="prefix for the knockout library read file path")

    parser.add_argument("--input_replicates", type=str, help="Input replicate indices")
    parser.add_argument("--control_replicates", type=str, help="Control replicate indices")
    parser.add_argument("--ko_replicates", type=str, help="KO replicate indices")

    parser.add_argument("--input_pairs", type=str, help="Input read pair indices")
    parser.add_argument("--control_pairs", type=str, help="Control read pair indices")
    parser.add_argument("--ko_pairs", type=str, help="KO read pair indices")

    parser.add_argument("--input_suffix", type=str, help="Input library read file suffix")
    parser.add_argument("--control_suffix", type=str, help="Control library read file suffix")
    parser.add_argument("--ko_suffix", type=str, help="KO library read file suffix")

    parser.add_argument("-u", "--umi", type=str, help="Library UMI indices", default="")

    parser.add_argument("--input_flag", action='store_false', help="Will not run the input pipeline")
    parser.add_argument("--control_flag", action='store_false', help="Will not run the control pipeline")
    parser.add_argument("--dedup_flag", action='store_false', help="Will not run the dedup pipeline")
    parser.add_argument("--align_flag", action='store_false', help="Will not run the alignment pipeline")
    parser.add_argument("--filter_flag", action='store_false', help="Will not run the filter pipeline")
    parser.add_argument("--peak_flag", action='store_false', help="Will not run the peak calling pipeline")
    parser.add_argument("--bigwig_flag", action='store_false', help="Will not run the bigwig creation pipeline")
    parser.add_argument("--cradle_flag", action='store_false', help="Will not run the cradle peak calling pipeline")

    args = parser.parse_args()

    input_lp = args.input_prefix
    control_lp = args.control_prefix
    ko_lp = args.ko_prefix

    input_lr = args.input_replicates
    control_lr = args.control_replicates
    ko_lr = args.ko_replicates

    input_lrp = args.input_pairs
    control_lrp = args.control_pairs
    ko_lrp = args.ko_pairs

    input_ls = args.input_suffix
    control_ls = args.control_suffix
    ko_ls = args.ko_suffix

    library_umis = args.umi
    if library_umis:
        all_umis = library_umis.strip().split()
        if len(all_umis)==3:
            input_library_umi_index = all_umis[0]
            control_library_umi_index = all_umis[1]
            ko_library_umi_index = all_umis[2]
        elif len(all_umis)==1:
            input_library_umi_index = all_umis[0]
            control_library_umi_index = all_umis[0]
            ko_library_umi_index = all_umis[0]
        else:
            raise ValueError("There should be three UMI indices for three libraries")

    # TODO: Reference genome and roi argument
    reference_genome = "/data5/deepro/genomes/hg38/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
    reference_genome_twobit = "/data5/deepro/genomes/hg38/hg38.2bit"
    roi_file = "/data5/deepro/starrseq/computational_pipeline/data/master.sorted.bed"
    
    finput=args.input_flag
    fcontrol=args.control_flag
    fdedup=args.dedup_flag
    falign=args.align_flag
    ffilter=args.filter_flag
    fpeak=args.peak_flag
    fbigwig=args.bigwig_flag
    fcpeak=args.cradle_flag

    main()
