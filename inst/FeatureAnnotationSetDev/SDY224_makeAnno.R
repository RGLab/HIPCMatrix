# Evan Henrich
# Sep 2018

library(ImmuneSpaceR)
con <- CreateConnection("SDY224", onTest = T)
gef <- con$getDataset("gene_expression_files")
gef <- gef[ gef$type == "PBMC", ]

library(GEOquery)
map <- lapply(gef$geo_accession, function(x){
  tmp <- getGEO(x)
  nm <- tmp@header$title
  nm <- gsub(" \\[(PBMC|Bcell)\\]", "", nm)
  return(c(x,nm))
})
map <- data.frame(t(data.frame(map)))
colnames(map) <- c("gsm","id")

# Get raw data from series matrix assuming named 'non-normalized.txt.gz' or 'raw.corrected.txt.gz'
baseDir <- paste0("/home/ehenrich/R/Gen_Test")
dir.create(baseDir)
gsm <- getGEO(gef[1, geo_accession])
gseList <- gsm@header$series_id
gseList <- gseList[ grepl("GSE45735", gseList)]
fls <- getGEOSuppFiles(gseList, makeDirectory = FALSE, baseDir = baseDir)
fls <- rownames(fls)
nms <- paste0("T", sub(".*T([0-9]+).*", "\\1", fls))
ems <- lapply(seq(1:length(fls)), function(x){
  flPath <- fls[[x]]
  nm <- nms[[x]]
  GEOquery::gunzip(flPath, overwrite = TRUE)
  flPath <- gsub(".gz", "", flPath)
  exprs <- data.table::fread(flPath, header = T)
  tmp <- colnames(exprs)[ grep("Gene", colnames(exprs), invert = T)]
  tmp <- paste0(nm, "_Day", tmp)
  colnames(exprs) <- c("Gene", tmp)
  colnames(exprs) <- as.character(map$gsm[ match(colnames(exprs), map$id)])
  colnames(exprs)[[1]] <- "feature_id"
  return(exprs)
})
exprs <- Reduce(function(x,y){ merge(x,y, by="feature_id", all=T)}, ems)
exprs <- exprs[ complete.cases(exprs)]

res <- data.frame(Probe_ID = exprs$feature_id,
                  Gene_Symbol = exprs$feature_id,
                  stringsAsFactors = F)
res <- res[ !is.na(res$Gene_Symbol), ]
res <- res[ !duplicated(res$Probe_ID), ]

# Write out
write.table(res,
            file = "FeatureAnnotationSetDev/SDY224_CustomAnno.tsv",
            quote = FALSE,
            sep = "\t",
            row.names = FALSE)

