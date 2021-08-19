# Create FAS called "HumanHT-12_v4" with note "made from Illumina manifest 2018".

manifest <- data.table::fread("FeatureAnnotationSetDev//HumanHT-12_V4_0_R2_15002873_B.txt",
                              skip = 8, fill = T)
manifest <- manifest[ !is.na(manifest$Entrez_Gene_ID), ]
manifest <- manifest[ , colnames(manifest) %in% c("Symbol", "Probe_Id"), with=F]
colnames(manifest) <- c("Gene_Symbol", "Probe_Id")
write.table(manifest,
            file = "FeatureAnnotationSetDev/IlluminaHumanHT-12_v4_2018.tsv",
            quote = FALSE,
            sep = "\t",
            row.names = FALSE)
