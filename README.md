# HIPCMatrix

Utilities for processing and analyzing [HIPC](https://www.immuneprofiling.org/hipc/page/show) gene expression data in [ImmuneSpace](https://www.immunespace.org/). 

## Processing directory structure: 

```
analysis_dir
├── analysis
│   └── exprs_matrices
│       ├── matrix_name.tsv
│       ├── matrix_name.tsv.raw
│       ├── matrix_name.tsv.summary
│       └── matrix_name.tsv.summary.orig
└── rawdata
    └── gene_expression
        ├── create-matrix
        │   └── matrix_name
        │       └── create-matrix-vars.tsv
        ├── RAW_FILES_FROM_IMMPORT
        └── supp_files
            └── sampletype_armaccession
                └── SDY1630_raw_expression.txt
```

* Raw files from ImmPort get saved into `analysis_dir/rawdata`
* Raw downloads from GEO or intermediate files from ImmPort are saved into
  `analysis_dir/rawdata/supp_files/sampletype_armaccession`. 
* Debug files get saved into `create-matrix/matrix_name`

## Full processing workflow

1. prep files
    1. pull files from correct location (GEO or supplementary files from ImmPort)
    1. If needed, clean up to be consistent per platform
        1. illumina: probe intensity and detection p-values for all samples 
           saved to one tsv file
        1. affymetrix: CEL files
        1. Stanford Genomics (two-color array): background-corrected expression 
           saved to one tsv file
        1. RNA-seq: raw counts saved to one tsv file
    1. Return path to cleaned file
1. prepare "raw" matrix: background-corrected expression for microarray; raw
   counts for rna-seq. Samples labeled using ImmPort biosample accession. 
    1. Load background-corrected or raw counts matrix from input file returned 
       from previous step
          1. illumina: Derive background-corrected expression using limma
          1. affymetrix: Derive background-corrected expression using RMA 
             (Result will be in log2 scale)
          1. two-color-array: Background correction already performed during 
             prep step, so only need to read input file
          1. RNA-seq: No transformation needed. Read raw counts into one table
    1. Map sample id to ImmPort biosample accession and subset to only include
       selected biosamples. 
    1. Return data.table with one row for feature, one column per biosample 
       accession. 
1. Normalize matrix: Normalize across samples to remove batch effects. Normalized
   matrix will be in log2 scale. 
    1. microarray: normalize across samples using quantile normalization. 
       Illumina data must first be log2 transformed. Affy and two color array 
       data already should be in log2 space based on raw matrix preparation
    1. rna-seq: DESeq2
1. Summarize matrix by gene symbol. This is the same process for all platforms
1. Write raw, normalized, and summarized matrices to disk as tsv
1. Write out log and debug files
    
    
    
