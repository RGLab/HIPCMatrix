# HIPCMatrix
<!-- badges: start -->
[![R build status](https://github.com/RGLab/HIPCMatrix/workflows/R-CMD-check/badge.svg)](https://github.com/RGLab/HIPCMatrix/actions)
[![Codecov test coverage](https://codecov.io/gh/RGLab/HIPCMatrix/branch/main/graph/badge.svg)](https://codecov.io/gh/RGLab/HIPCMatrix?branch=main)
<!-- badges: end -->

Utilities for processing and analyzing [HIPC](https://www.immuneprofiling.org/hipc/page/show) gene expression data in [ImmuneSpace](https://www.immunespace.org/). 

# Processing Gene Expression Data

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
    
# Gene Expression Analysis

Methods for downstream analysis used in ImmuneSpace are also maintained in this 
package. They are accessible through the `HMX` object, which is an R6 object 
that extends the `ISConn` object defined in [`ImmuneSpaceR`](https://github.com/RGLab/ImmuneSpaceR). 

To create an `HMX` object: 

```
con <- HMX$new("SDY269")
```

You can then perform analyses on the matrices available through the connection 
object. For example, to run differential expression analysis: 

```
de_results <- con$runGEAnalysis()
```

# Feature annotation

This package provides tools for ensuring that the Gene Symbols used throughout
the ImmuneSpace portal are consistent and up to date.  The `UpdateAnno` package 
previously managed gene annotation in ImmuneSpace. 
There are a number of places where Gene Symbols play an important role:
- Gene Expression Explorer - Pulls Gene Symbols from the FeatureAnnotation query
- Gene Set Enrichment Analysis - Pulls Gene Symbols from the static gene set 
module data objects in this package, e.g. chaussabel or emory.
- Differential Gene Expression Analysis - Uses the GeneExpressionAnalysisResults 
query to display differentially expressed genes.

## Custom Feature Annotation

When loading a new matrix onto ImmuneSpace, it must be associated with a 
feature annotation set, saved in the Feature Annotation Set table. If the 
annotation for the matrix is not already loaded from a previous study which 
used the same platform, you must create a new one. Each feature annotation 
table and the source code to create it is saved in `HIPCMatrix/inst/FeatureAnnotationSetDev`. 

Please refer to the 
[Notion Documentation](https://www.notion.so/rglab/Load-a-new-gene-expression-matrix-385ce687af594d369554e406864c12ad) 
for more details on how to create and upload a FAS. 

## Updating annotation throughout ImmuneSpace 

The function `create_gene_alias_map()` pulls the most recent annotation from 
the latest HUGO Gene Nomenclature Committee dataset. 
This resource is used instead of the NCBI database via biomaRt or the equivalent 
`org.Hs.eg.db` as those data sources contained mappings that were deemed incorrect 
(such as "ACTB" > "POTEF") during the ImmuneSignatures 2 project.

This function handles 2 edge cases: 
1. An alias maps to itself as a symbol as well as other symbols.  In this case, 
we have selected to use the self-mapping and removed any other mappings.  
Our rationale is that the other-mapping symbols are historic artifacts and no 
longer accurate.
2. An alias maps to multiple symbols that do not include itself.  In this case, 
we drop these aliases from the matrix because we do not have a good way to know 
which symbol is the most accurate. 

This mapping table, as well as gene sets with the most current mapping, are 
saved in HIPCMatrix package data to be easily accessible from many locations. 

The `HMX` object includes methods for updating annotation in ImmuneSpace, 
using the gene symbol mapping table installed with the package: 

To update feature annotation sets associated with matrices associated with a 
connection: 
```
con <- HMX$new("")
con$updateFAS()
```

To update summary.tsv expressionsets associated with a connection: 
```
con$updateEMs()
```



To update these mappings and apply it throughout ImmuneSpace: 

1. Source `HIPCMatrix/data-raw/update_hgnc_mapping.R`. 
    This will update the HGNC mapping table which is part of the HIPCMatrix 
    package data, then use this new mapping to update the gene signatures 
    which are also part of the package data. 
1. Bump package version
1. Commit the changes, and push to github
1. Install the update on the server: 
    ```
    ssh rsT
    s
    R
    devtools::install("RGLab/HIPCMatrix@dev")
    ```
1. On the server, go to data integration module at Studies level: 
    * TEST: https://test.immunespace.org/dataintegration/Studies/begin.view?
    * PROD: https://www.immunespace.org/dataintegration/Studies/begin.view?

1. Run Update Anno ETL. This will: 
    1. Update Feature Annotation Sets with current annotation
    1. Update .summary.tsv expressionsets with the current annotation
    1. Re-run GSEA for all studies where it is turned on. 
