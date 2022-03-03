#!/share/apps/R/R-4.0.2/bin/Rscript



setwd("/home/csexton/compute/NMF_differentiation/RcppML_nmf")

# devtools::install_github("zdebruine/RcppSparse") # compile dev version
# devtools::install_github("zdebruine/RcppML") # compile dev version
# install.packages("RcppML")

library(Matrix)
library(RcppML)
library(tidyverse)
library(pheatmap)

con_mat <- read_tsv("matrix_connections.tsv")

x <- as.matrix(con_mat[,-1])
rownames(x) <- con_mat$region
matCSC <- as(x, "dgCMatrix")

num_labels = 3
num_runs = 100

model <- nmf(matCSC, k = num_labels, maxit = num_runs, tol = 1e-100)

# saveRDS(model,paste0("model_",num_labels,".rds"))
#colnames(model$h) <- colnames(matCSC)
#rownames(model$h) <- 1:num_labels

png(paste0("h_",num_labels,".png"),width = 1300, height = 1000, res = 200)
pheatmap(t(model$h), scale = "row", border_color = NA, cluster_rows = F,
         cluster_cols = T, fontsize = 14,
         main = paste0("H matrix k=",num_labels),
         breaks = seq(-1, 1, length.out = 100))
# main = "Segmentation of H1 and naive hESCs",
# breaks = seq(-1, 1, length.out = 100)
dev.off()

#colnames(model$w) <- 1:num_labels
#rownames(model$w) <- rownames(matCSC)


z <- as.data.frame(ifelse(grepl("endo",rownames(matCSC)),"endo","h1"))
rownames(z) <- rownames(matCSC)
colnames(z) <- "cell_type"

# small <- model$w[c(1:50,20000:20050),]

png(paste0("w_",num_labels,".png"),width = 1300, height = 1000, res = 200)
pheatmap(model$w, scale = "row", border_color = NA, cluster_rows = T,
         cluster_cols = F, fontsize = 14, annotation_row = z,show_rownames = F,
         main = paste0("W matrix k=",num_labels))
dev.off()






