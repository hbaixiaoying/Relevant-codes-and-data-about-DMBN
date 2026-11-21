
### HW_A differential multilayer biomolecular network framework for decoding 
### multi-omics perturbations in lung squamous cell carcinoma
rm(list = ls())
library(tidyverse)
library(psych)
library(data.table)

setwd(".../DMBN_code")

### ---- Load raw data ----
mirna   <- read.csv("Primary_data/mirna_mature_tpm_log2_row_norm.csv")
rna     <- read.csv("Primary_data/rna_fpkm_log2_NArm_row_norm.csv")
protein <- read.csv("Primary_data/global_proteome_ratio_norm_NArm_gene_level.csv")

### ---- Drop annotation columns, keep one ID column ----
mirna   <- mirna[,  -c(1, 2, 3, 5, 6)]
rna     <- rna[,    -c(1, 2, 3, 5)]
protein <- protein[, -c(1:17)]

### miRNA: drop features with >70% NA (ID column excluded from denominator);
### RNA/protein were pre-filtered upstream, so only impute remaining NAs.
mirna <- mirna[rowSums(is.na(mirna)) / (ncol(mirna) - 1) <= 0.7, ]
mirna[is.na(mirna)]     <- 0
rna[is.na(rna)]         <- 0
protein[is.na(protein)] <- 0

### ---- Collapse duplicate IDs by mean, move ID to rownames ----
mirna   <- aggregate(. ~ Name,       data = mirna,   mean) %>% column_to_rownames("Name")
protein <- aggregate(. ~ GeneSymbol, data = protein, mean) %>% column_to_rownames("GeneSymbol")
rna     <- rna %>% column_to_rownames(names(rna)[1])

### ---- Match tumor/normal samples across all omics ----
samples    <- Reduce(intersect, list(colnames(rna), colnames(mirna), colnames(protein)))

### Keep only sample IDs with both tumor and paired normal samples.
normal_ids <- grep("\\.N$", samples, value = TRUE)
tumor_ids  <- sub("\\.N$", "", normal_ids)
paired     <- tumor_ids %in% samples
tumor_ids  <- tumor_ids[paired]
normal_ids <- normal_ids[paired]

### ---- Split each omics matrix into tumor and normal ----
tumor_rna     <- rna[,     tumor_ids,  drop = FALSE]
tumor_mirna   <- mirna[,   tumor_ids,  drop = FALSE]
tumor_protein <- protein[, tumor_ids,  drop = FALSE]

normal_rna     <- rna[,     normal_ids, drop = FALSE]
normal_mirna   <- mirna[,   normal_ids, drop = FALSE]
normal_protein <- protein[, normal_ids, drop = FALSE]

### ---- Write tumor and normal matrices ----
#write.csv(tumor_rna,     "tumor_rna.csv",     row.names = TRUE)
#write.csv(tumor_mirna,   "tumor_mirna.csv",   row.names = TRUE)
#write.csv(tumor_protein, "tumor_protein.csv", row.names = TRUE)
#write.csv(normal_rna,     "normal_rna.csv",     row.names = TRUE)
#write.csv(normal_mirna,   "normal_mirna.csv",   row.names = TRUE)
#write.csv(normal_protein, "normal_protein.csv", row.names = TRUE)
