# StarrSeq peak calling pipeline

## Conda environment for starrseq data exploration and analysis
```bash
foo@bar:~$ conda create -n starrseq -c conda-forge -c bioconda python=3.9 matplotlib jupyter ipykernel bwa samtools picard bedtools pandas openpyxl scikit-learn multiprocess bokeh
```

## Conda environment for starrpeaker

```bash
foo@bar:~$ conda create -n starrpeaker -c bioconda python=2.7 pybedtools
foo@bar:~$ conda activate starrpeaker
foo@bar:~$ pip install git+https://github.com/gersteinlab/starrpeaker
```
**STARRPEAKER PROBLEM: https://stackoverflow.com/questions/43147475/when-i-use-for-loop-indexerror-list-index-out-of-range**

## Conda environment for cradle

```bash
foo@bar:~$ conda create -n cradle python=3.7
foo@bar:~$ conda activate cradle
foo@bar:~$ pip install synapseclient
foo@bar:~$ pip install git+https://github.com/deeprob/cradle
```

## Naming convention of libraries
${LIBRARYPREFIX}_${REPLICATEINDEX}_${READPAIRID}_${LIBRARYSUFFIX}

## Align existing reads to the reference genome GRCh38 and create BAM files

- Step 1: Download the reference genome
    ```bash
    foo@bar:~$ wget https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz
    foo@bar:~$ gunzip GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz
    ```
- Step 2: Create reference index using bwa with -a bwtsw option for human genome
    ```bash
    foo@bar:~$ bwa index -a bwtsw GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
    ```
- Step 3: Align paired end reads to the reference genome using bwa mem and convert them to bam format

    **Read files come from both input and output libraries (STARRSeq naming conventions). Each library will have three replicates and two read files per replicate. Can use the align_fastq_to_ref.sh script for aligning libraries.**

## Filter reads
- Step 1: Use picard to get rid of duplicates in the aligned BAM files if there are no UMIs, else use STARRDUST before alignment -  estimated time: 3 hours

- Step 2: Use samtools to filter unmapped, secondary alignments, mapping quality score less than 40 - estimated time: few minutes

- Step 3: Pool the replicates in a library for calling peaks - estimated time: few seconds

**NB: Avoid piping while using samtools.**

## Peak calling using starrpeaker

- Step 1: Download and unzip covariate files given in starrpeaker github - estimated time: ~30 mins
    Downloaded from the link: https://helkar.synology.me:5001/fsdownload/EKNVQrVHM/covariates
    **NB: Could not download precomputed covariates file using wget. Had to download through web browser and scp to cluster.** 
- Step 2: Download the chromosome size file - estimated time: few seconds
    ```bash
    foo@bar:~$ wget https://raw.githubusercontent.com/gersteinlab/starrpeaker/master/data/GRCh38.chrom.sizes.simple.sorted
    ```
- Step 3: Download the blacklisted regions file - estimated time: few seconds
    ```bash
    foo@bar:~$ wget https://raw.githubusercontent.com/gersteinlab/starrpeaker/master/data/ENCODE_blacklist_GRCh38_ENCFF419RSJ_merged.bed
    ```
- Step 4: Peak calling - estimated time: ~2.5 hours
