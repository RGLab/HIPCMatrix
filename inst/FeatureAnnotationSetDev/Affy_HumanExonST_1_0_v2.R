# SDY1291 - custom cdf from GPL10647 on GEO

# Following is script used to create package, stored as a tar.gz in UpdateAnno
# pkgName: huex10stv2cdf
# AnnoNm: HuEx-1_0-st-v2
# BiocInstaller::biocLite("makecdfenv")
library(makecdfenv)
setwd("/")
makecdfenv::make.cdf.package(filename = "GPL10647.cdf",
                             packagename = "huex10stv2cdf",
                             cdf.path = "home/ehenrich/Downloads/",
                             package.path = "/home/ehenrich/R/UpdateAnno/FeatureAnnotationSetDev/",
                             species = "Homo_sapiens",
                             compress = TRUE)

# create tarball in same directory
tar(tarfile = "/home/ehenrich/R/UpdateAnno/FeatureAnnotationSetDev/huex10stv2cdf.tar.gz",
    files = "/home/ehenrich/R/UpdateAnno/FeatureAnnotationSetDev/huex10stv2cdf/",
    compression = "gzip")

# test on Cel file from SDY1291 and also creating custom annotation:
# untar("/home/ehenrich/R/UpdateAnno/FeatureAnnotationSetDev/huex10stv2cdf.tar.gz")
# R CMD INSTALL /FeatureAnnotationSetDev/huex10stv2cdf
inputFiles <- "/home/ehenrich/Downloads/GSM563106.CEL"
setwd("/") # b/c filepaths are absolute and justRMA prepends wd
eset <- affy::justRMA(filenames = inputFiles, normalize = FALSE)
em <- Biobase::exprs(eset)

em <- data.frame(em, stringsAsFactors = F)
affyids <- rownames(em)
entrezids <- gsub("_at", "", affyids)
library(org.Hs.eg.db)
library(annotate)
library(data.table)
std <- suppressMessages(data.table(AnnotationDbi::select(org.Hs.eg.db,
                                                         keys = AnnotationDbi::keys(org.Hs.eg.db,
                                                                                    keytype = "SYMBOL"),
                                                         columns = c("ENTREZID", "SYMBOL"),
                                                         keytype = "SYMBOL")))
gs <- std$SYMBOL[ match(entrezids, std$ENTREZID)]

res <- data.frame(Probe_ID = affyids,
                  Gene_Symbol = gs,
                  stringsAsFactors = F)
res <- res[ !is.na(res$Gene_Symbol), ]
res <- res[ !duplicated(res$Probe_ID), ]

# Write out
write.table(res,
            file = "/home/ehenrich/R/UpdateAnno/FeatureAnnotationSetDev/Affy_HumanExonST_1_0_v2.tsv",
            quote = FALSE,
            sep = "\t",
            row.names = FALSE)
