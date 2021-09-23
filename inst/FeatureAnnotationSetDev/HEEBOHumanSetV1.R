# E.henrich
# Jan 2019
# HEEBO Human Set V1 is a two-color array created by the Stanford Functional Genomics Facility.
# The probe to symbol mapping is available as an .xls at
# https://microarray.org/data/download/HEEBO_Human_Set_v1.00.xls
# The xls was downloaded and then edited in libreoffice to create a tsv that is saved
# here as featureAnnotationSetDev/HEEBOHumanSetV1_unformatted.tsv

tmp <- data.table::fread("FeatureAnnotationSetDev/HEEBOHumanSetV1_unformatted.tsv")

# Remove unmapped
res <- data.frame(Probe_ID = tmp$Oligo_ID,
                  Gene_Symbol = tmp$`Symbol v12`,
                  stringsAsFactors = F)
res <- res[ !is.na(res$Gene_Symbol), ]

# Check all are unique
res <- res[ !duplicated(res), ]

# Remove unnamed genes
res <- res[ !grepl("(LOC\\d+|^\\.$)", res$Gene_Symbol), ]

# Remove misnamed genes (converted by xls to date - found this way in orig file)
bad <- grepl("\\d{1,2}-[[:alpha:]]{3}", res$Gene_Symbol)
res <- res[ !bad, ]

# Write out
write.table(res,
            file = "FeatureAnnotationSetDev/HEEBOHumanSetV1.tsv",
            quote = FALSE,
            sep = "\t",
            row.names = FALSE)
