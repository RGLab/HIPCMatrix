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
  1. 
