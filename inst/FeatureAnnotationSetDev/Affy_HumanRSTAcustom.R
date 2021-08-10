# SDY1328 - custom cdf from GPL15048 on GEO

# Following is script used to create package, stored as a tar.gz in UpdateAnno
# Would need to re-download the necessary .CDF and .CEL files from GPL15048 and GSM1607320

# pkgName: hursta2a520709cdf
# AnnoNm: HuRSTA-2a520709
# BiocInstaller::biocLite("makecdfenv")
library(makecdfenv)
setwd("/")
makecdfenv::make.cdf.package(filename = "GPL15048_HuRSTA_2a520709.CDF",
                             packagename = "hursta2a520709cdf",
                             cdf.path = "/home/ehenrich/R/UpdateAnno/FeatureAnnotationSetDev/",
                             package.path = "/home/ehenrich/R/UpdateAnno/FeatureAnnotationSetDev/",
                             species = "Homo_sapiens",
                             compress = TRUE)

# create tarball in same directory
tar(tarfile = "/home/ehenrich/R/UpdateAnno/FeatureAnnotationSetDev/hursta2a520709cdf.tar.gz",
    files = "/home/ehenrich/R/UpdateAnno/FeatureAnnotationSetDev/hursta2a520709cdf/",
    compression = "gzip")

# test on Cel file from SDY1291 and also creating custom annotation:
# untar("/home/ehenrich/R/UpdateAnno/FeatureAnnotationSetDev/hursta2a520709cdf.tar.gz")
# R CMD INSTALL /FeatureAnnotationSetDev/hursta2a520709cdf
inputFiles <- "/home/ehenrich/R/UpdateAnno/FeatureAnnotationSetDev/GSM1607320__SDY1328smpl.CEL"
setwd("/") # b/c filepaths are absolute and justRMA prepends wd
eset <- affy::justRMA(filenames = inputFiles, normalize = FALSE)
em <- Biobase::exprs(eset)
em <- data.frame(em, stringsAsFactors = F)
affyids <- rownames(em)

# get mappings from GPL15048
library(GEOquery)
tmp <- getGEO("GPL15048")
mp <- tmp@dataTable@table
library(org.Hs.eg.db)
library(annotate)
library(data.table)
std <- suppressMessages(data.table(AnnotationDbi::select(org.Hs.eg.db,
                                                         keys = AnnotationDbi::keys(org.Hs.eg.db,
                                                                                    keytype = "SYMBOL"),
                                                         columns = c("ENTREZID", "SYMBOL"),
                                                         keytype = "SYMBOL")))
gs <- std$SYMBOL[ match(mp$EntrezGeneID, std$ENTREZID)]

res <- data.frame(Probe_ID = affyids,
                  Gene_Symbol = gs,
                  stringsAsFactors = F)
res <- res[ !is.na(res$Gene_Symbol), ]
res <- res[ !duplicated(res$Probe_ID), ]

# Write out
write.table(res,
            file = "/home/ehenrich/R/UpdateAnno/FeatureAnnotationSetDev/Affy_HumanRSTAcustom.tsv",
            quote = FALSE,
            sep = "\t",
            row.names = FALSE)
