# STARRSuite
A set of tools to analyze data produced by STARRSeq libraries

## Creating conda environments 

### Base conda environment for data exploration and analysis

```bash
foo@bar:~$ conda create -n starrseq -c conda-forge -c bioconda -c anaconda -c bjrn python=3.9 matplotlib jupyter ipykernel bwa samtools picard pandas openpyxl scikit-learn multiprocess bokeh requests pybedtools pysam bioinfokit matplotlib-venn adjusttext
```

### Conda environment for starrpeaker peak calling

```bash
foo@bar:~$ conda create -n starrpeaker -c bioconda python=2.7 pybedtools
foo@bar:~$ conda activate starrpeaker
foo@bar:~$ pip install git+https://github.com/deeprob/starrpeaker
```

### Conda environment for cradle peak calling

```bash
foo@bar:~$ conda create -n cradle python=3.7
foo@bar:~$ conda activate cradle
foo@bar:~$ pip install synapseclient deeptools
foo@bar:~$ pip install git+https://github.com/deeprob/cradle
```

### Conda environment for macs2 peak calling

```bash
foo@bar:~$ conda create -n macs2 -c conda-forge -c bioconda macs2 -y
```

## Downloading datasets required by the tools

### Reference genome GRCh38

- Step 1: Download the reference genome
    ```bash
    foo@bar:~$ wget https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz
    foo@bar:~$ gunzip GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz
    ```
- Step 2: Create reference index using bwa with -a bwtsw option for human genome
    ```bash
    foo@bar:~$ bwa index -a bwtsw GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta
    ```

### STARRPeaker data

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

### CRADLE data

TODO

### 
