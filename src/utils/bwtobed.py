import os
import pyBigWig


bw_file = "/data5/deepro/encode/hek293/atacseq/ENCFF668ATM.bigWig"
bed_file = os.path.splitext(bw_file)[0] + ".bed"
bw = pyBigWig.open(bw_file)  
of = open(bed_file, "w")

for chrom, len in bw.chroms().items():
    intervals = bw.intervals(chrom)
    for interval in intervals:
        of.write(f"{chrom}\t{interval[0]}\t{interval[1]}\t{interval[2]}\n")


bw.close()
of.close()