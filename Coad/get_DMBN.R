
### HW_A differential multilayer biomolecular network framework for decoding 
### multi-omics perturbations in lung squamous cell carcinoma
rm(list=ls())
library(tidyverse)
library(psych)
library(data.table)
library(igraph)
library(Matrix)
setwd(".../DMBN_code")

mirna_n<-read.csv("Output/normal_mirna.csv",row.names = 1)
mirna_t<-read.csv("Output/tumor_mirna.csv",row.names = 1)

rna_n<-read.csv("Output/normal_rna.csv",row.names = 1)
rna_t<-read.csv("Output/tumor_rna.csv",row.names = 1)

protein_n<-read.csv("Output/normal_protein.csv",row.names = 1)
protein_t<-read.csv("Output/tumor_protein.csv",row.names = 1)

mirna_n<-as.matrix(t(mirna_n))
mirna_t<-as.matrix(t(mirna_t))

rna_n<-as.matrix(t(rna_n))
rna_t<-as.matrix(t(rna_t))

protein_n<-as.matrix(t(protein_n))
protein_t<-as.matrix(t(protein_t))

mirna_n_cor <- cor(mirna_n, use = "all.obs",method="spearman")
mirna_t_cor <- cor(mirna_t, use = "all.obs",method="spearman")

rna_n_cor <- cor(rna_n, use = "all.obs",method="spearman")
rna_t_cor <- cor(rna_t, use = "all.obs",method="spearman")

protein_n_cor <- cor(protein_n, use = "all.obs",method="spearman")
protein_t_cor <- cor(protein_t, use = "all.obs",method="spearman")

### Cutoffs are the |rho| value at the top 0.05% (5e-4) of all unique
### edge correlations (|rho| ranking), computed separately per omics layer and group.

mirna_n_adj <- abs(mirna_n_cor) >= 0.6945481
mirna_t_adj <- abs(mirna_t_cor) >= 0.7037044
diff_mirna_adj<-abs(mirna_n_adj - mirna_t_adj)

rna_n_adj <- abs(rna_n_cor) >= 0.8256872
rna_t_adj <- abs(rna_t_cor) >= 0.7266718
diff_rna_adj<-abs(rna_n_adj - rna_t_adj)

protein_n_adj <- abs(protein_n_cor) >= 0.6894900
protein_t_adj <- abs(protein_t_cor) >= 0.6887888
diff_protein_adj<-abs(protein_n_adj - protein_t_adj)


diff_mirna_adj <- diff_mirna_adj[
  rowSums(diff_mirna_adj != 0) > 0,
  colSums(diff_mirna_adj != 0) > 0,
  drop = FALSE
]

diff_rna_adj <- diff_rna_adj[
  rowSums(diff_rna_adj != 0) > 0,
  colSums(diff_rna_adj != 0) > 0,
  drop = FALSE
]

diff_protein_adj <- diff_protein_adj[
  rowSums(diff_protein_adj != 0) > 0,
  colSums(diff_protein_adj != 0) > 0,
  drop = FALSE
]


### ---- 1. Convert protein node names to numeric IDs ----
protein_id_map <- data.frame(
  Original = rownames(diff_protein_adj),
  New      = as.character(seq_len(nrow(diff_protein_adj))),
  stringsAsFactors = FALSE
)

dimnames(diff_protein_adj) <- list(protein_id_map$New, protein_id_map$New)

rna_nodes <- rownames(diff_rna_adj)
shared_ids <- intersect(rna_nodes, protein_id_map$Original)
protein_id_map_shared <- protein_id_map[protein_id_map$Original %in% shared_ids, ]


mirna_rna_edge<-read.csv("Primary_data/mirna_rna_interlayer_edge_targetscan.csv")
rna_protein_edge<-protein_id_map_shared
protein_rna_edge<-protein_id_map_shared[, c(2, 1)]

### ---- Node order: miRNA, RNA, protein ----
mirna_nodes   <- rownames(diff_mirna_adj)
rna_nodes     <- rownames(diff_rna_adj)
protein_nodes <- rownames(diff_protein_adj)   # already numeric IDs
all_nodes     <- c(mirna_nodes, rna_nodes, protein_nodes)

n <- length(all_nodes)
DMBN <- matrix(0, n, n, dimnames = list(all_nodes, all_nodes))

### ---- Diagonal blocks ----
DMBN[mirna_nodes,   mirna_nodes]   <- diff_mirna_adj
DMBN[rna_nodes,     rna_nodes]     <- diff_rna_adj
DMBN[protein_nodes, protein_nodes] <- diff_protein_adj


### ---- miRNA -> RNA block ----
idx <- cbind(mirna_rna_edge$mirna, mirna_rna_edge$rna)
DMBN[idx] <- 1

### ---- RNA -> protein block ----
idx <- cbind(protein_id_map_shared$Original, protein_id_map_shared$New)
DMBN[idx] <- 1

### ---- protein -> RNA block ----
idx <- cbind(protein_id_map_shared$New, protein_id_map_shared$Original)
DMBN[idx] <- 1

save(DMBN, mirna_nodes, rna_nodes, protein_nodes, all_nodes, protein_id_map,
     file = "Output/get_DMBN_objects.RData")
