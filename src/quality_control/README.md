# STARRSeq library report
This package generates a bunch of quality control metrics for STARRSeq input and output libraries

## Read Quality Control Metrics

1. Number of Regions of Interest (ROIs)
2. Mean Size of ROIs
3. Total ROI size in bps
4. Number of Raw Reads (All replicates)
5. Number of Filtered Reads (All replicates)
6. Sequencing Read Length 
7. Type of sequencing (SE or PE)
8. Coverage (per replicate and merged)
9. Number of ROIs with at least N reads (per replicate and merged)
10. Percentage of ROIs with at least N reads (per replicate and merged)

## Library Quality Control Metrics

1. Input-input replicate correlation
2. Output-output replicate correlation
3. Input-output replicate correlation
4. Enhancer summit correlation between replicates
5. Fold change correlation between replicates
6. Reproducible peaks between replicates

## Methods transparency

1. Sequencer
2. Demultiplexing tool and commands
3. Raw read qc tool and report
4. Deduplication tool and command
5. Alignment tool and command
6. Filter tool and command
7. Peak calling tool and command
