# Evan Henrich
# September 2019
# Notes: SDY1256 / SDY1412 are RNASeq study using same ensembl build
# example accession: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3030156

library(data.table)
library(biomaRt)

# Generate mapping using build found on GSM 'data processing' notes: Ensembl GRCh38.86
# Get ensembl id and gene id
mart <- useEnsembl(biomart = "ensembl", version = 86, dataset = "hsapiens_gene_ensembl")
ens2gene <- getBM(attributes=c('ensembl_gene_id','hgnc_symbol'),
                  mart = mart)

# Remove unmapped and duplicated probe ids
res <- data.frame(Probe_ID = ens2gene$ensembl_gene_id,
                  Gene_Symbol = ens2gene$hgnc_symbol,
                  stringsAsFactors = F)
res <- res[ !is.na(res$Gene_Symbol) & res$Gene_Symbol != "", ]
res <- res[ !duplicated(res$Probe_ID), ]

# Write out
write.table(res,
            file = "FeatureAnnotationSetDev/Ensembl_86_GRCh38.tsv",
            quote = FALSE,
            sep = "\t",
            row.names = FALSE)



