
### HW_A differential multilayer biomolecular network framework for decoding 
### multi-omics perturbations in lung squamous cell carcinoma
rm(list = ls())
library(igraph)
library(Matrix)

load("Output/get_DMBN_objects.RData")

g <- graph_from_adjacency_matrix(DMBN, mode = "directed", weighted = NULL)

### ---- Directed degree centrality: total = in + out ----
deg <- degree(g, mode = "total", loops = FALSE, normalized = FALSE)

### ---- Directed betweenness centrality ----
betw <- betweenness(g, directed = TRUE, weights = NA, normalized = FALSE)

### ---- Align to all_nodes order ----
deg  <- deg[all_nodes]
betw <- betw[all_nodes]
### ---- L1 normalization -> P0 ----
normalize_l1 <- function(x) {
  if (sum(x) == 0) return(rep(1/length(x), length(x)))
  x / sum(x)
}

deg_norm  <- normalize_l1(deg)
betw_norm <- normalize_l1(betw)
P0 <- normalize_l1(deg_norm + betw_norm)


### ---- Slice DMBN into 9 blocks ----
n_m <- length(mirna_nodes)
n_r <- length(rna_nodes)
n_p <- length(protein_nodes)

miRNA_idx   <- 1:n_m
RNA_idx     <- (n_m + 1):(n_m + n_r)
Protein_idx <- (n_m + n_r + 1):(n_m + n_r + n_p)

network_blocks <- list(
  mm = DMBN[miRNA_idx,   miRNA_idx,   drop = FALSE],
  rr = DMBN[RNA_idx,     RNA_idx,     drop = FALSE],
  pp = DMBN[Protein_idx, Protein_idx, drop = FALSE],
  mr = DMBN[miRNA_idx,   RNA_idx,     drop = FALSE],
  mp = DMBN[miRNA_idx,   Protein_idx, drop = FALSE],
  rm = DMBN[RNA_idx,     miRNA_idx,   drop = FALSE],
  rp = DMBN[RNA_idx,     Protein_idx, drop = FALSE],
  pm = DMBN[Protein_idx, miRNA_idx,   drop = FALSE],
  pr = DMBN[Protein_idx, RNA_idx,     drop = FALSE]
)

### ---- Multilayer RWR transition matrix (from old script, unchanged) ----
multilayer_RWR_strict <- function(adj_blocks, lambda = 0.7) {
  n_mirna    <- nrow(adj_blocks$mm)
  n_rna      <- nrow(adj_blocks$rr)
  n_protein  <- nrow(adj_blocks$pp)
  total_nodes <- n_mirna + n_rna + n_protein
  
  all_nodes <- c(rownames(adj_blocks$mm),
                 rownames(adj_blocks$rr),
                 rownames(adj_blocks$pp))
  
  H <- Matrix(0, nrow = total_nodes, ncol = total_nodes, sparse = TRUE)
  rownames(H) <- all_nodes
  colnames(H) <- all_nodes
  
  miRNA_idx   <- 1:n_mirna
  RNA_idx     <- (n_mirna + 1):(n_mirna + n_rna)
  Protein_idx <- (n_mirna + n_rna + 1):total_nodes
  
  norm_vec <- function(x) {
    s <- sum(x)
    if (s == 0) return(x)
    x / s
  }
  
  ### miRNA layer
  for (i in seq_len(n_mirna)) {
    k <- 0
    if (sum(adj_blocks$mr[i, ]) > 0) k <- k + 1
    if (sum(adj_blocks$mp[i, ]) > 0) k <- k + 1
    if (k > 0) {
      H[miRNA_idx[i], miRNA_idx] <- (1 - lambda) * norm_vec(adj_blocks$mm[i, ])
    } else {
      H[miRNA_idx[i], miRNA_idx] <- norm_vec(adj_blocks$mm[i, ])
    }
    if (k > 0) {
      w <- lambda / k
      if (sum(adj_blocks$mr[i, ]) > 0)
        H[miRNA_idx[i], RNA_idx] <- w * norm_vec(adj_blocks$mr[i, ])
      if (sum(adj_blocks$mp[i, ]) > 0)
        H[miRNA_idx[i], Protein_idx] <- w * norm_vec(adj_blocks$mp[i, ])
    }
  }
  
  ### RNA layer
  for (i in seq_len(n_rna)) {
    k <- 0
    if (sum(adj_blocks$rm[i, ]) > 0) k <- k + 1
    if (sum(adj_blocks$rp[i, ]) > 0) k <- k + 1
    if (k > 0) {
      H[RNA_idx[i], RNA_idx] <- (1 - lambda) * norm_vec(adj_blocks$rr[i, ])
    } else {
      H[RNA_idx[i], RNA_idx] <- norm_vec(adj_blocks$rr[i, ])
    }
    if (k > 0) {
      w <- lambda / k
      if (sum(adj_blocks$rm[i, ]) > 0)
        H[RNA_idx[i], miRNA_idx] <- w * norm_vec(adj_blocks$rm[i, ])
      if (sum(adj_blocks$rp[i, ]) > 0)
        H[RNA_idx[i], Protein_idx] <- w * norm_vec(adj_blocks$rp[i, ])
    }
  }
  
  ### Protein layer
  for (i in seq_len(n_protein)) {
    k <- 0
    if (sum(adj_blocks$pm[i, ]) > 0) k <- k + 1
    if (sum(adj_blocks$pr[i, ]) > 0) k <- k + 1
    if (k > 0) {
      H[Protein_idx[i], Protein_idx] <- (1 - lambda) * norm_vec(adj_blocks$pp[i, ])
    } else {
      H[Protein_idx[i], Protein_idx] <- norm_vec(adj_blocks$pp[i, ])
    }
    if (k > 0) {
      w <- lambda / k
      if (sum(adj_blocks$pm[i, ]) > 0)
        H[Protein_idx[i], miRNA_idx] <- w * norm_vec(adj_blocks$pm[i, ])
      if (sum(adj_blocks$pr[i, ]) > 0)
        H[Protein_idx[i], RNA_idx] <- w * norm_vec(adj_blocks$pr[i, ])
    }
  }
  return(H)
}

### ---- PRWR-H iteration (from old script, unchanged) ----
rwr_iteration <- function(transition_matrix, P0, restart_prob = 0.5,
                          max_iter = 3000, tol = 1e-6) {
  n <- length(P0)
  p <- P0
  for (i in seq_len(max_iter)) {
    p_next <- (1 - restart_prob) * (t(transition_matrix) %*% p) + restart_prob * p
    p_next <- as.vector(p_next)
    if (sum(abs(p_next - p)) < tol) {
      cat("Converged at iteration", i, "\n")
      return(p_next)
    }
    p <- p_next
  }
  cat("Reached max iteration", max_iter, "\n")
  return(p)
}

### ---- Build transition matrix and run PRWR-H ----
transition_matrix <- multilayer_RWR_strict(network_blocks, lambda = 0.7)

### Run PRWR-H
result <- rwr_iteration(transition_matrix, P0, restart_prob = 0.5)

### ---- Output ----
results_df <- data.frame(
  Node  = rownames(transition_matrix),
  Score = result
)

results_df$Layer <- ifelse(results_df$Node %in% mirna_nodes,   "miRNA",
                           ifelse(results_df$Node %in% rna_nodes,     "RNA",
                                  ifelse(results_df$Node %in% protein_nodes, "protein", NA)))

mirna_rwr   <- results_df[results_df$Layer == "miRNA",   c("Node", "Score")]
rna_rwr     <- results_df[results_df$Layer == "RNA",     c("Node", "Score")]
protein_rwr <- results_df[results_df$Layer == "protein", c("Node", "Score")]

### ---- Sort each layer by Score (descending) ----
mirna_rwr2   <- mirna_rwr[order(-mirna_rwr$Score), ]
rna_rwr2     <- rna_rwr[order(-rna_rwr$Score), ]
protein_rwr2 <- protein_rwr[order(-protein_rwr$Score), ]

### ---- Add original protein symbol back via protein_id_map ----
protein_rwr2$Original <- protein_id_map$Original[match(protein_rwr2$Node, protein_id_map$New)]

save(
  results_df,                              
  mirna_nodes, rna_nodes, protein_nodes, all_nodes,  
  deg, betw, eigen,                       
  protein_id_map,                          
  DMBN,                                   
  file = "Output/get_PRWR-H_objects.RData"
)

