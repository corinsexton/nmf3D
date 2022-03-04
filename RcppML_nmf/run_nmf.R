#!/share/apps/R/R-4.0.2/bin/Rscript



# setwd("/home/csexton/compute/NMF_differentiation/RcppML_nmf")
setwd("/Users/coripenrod/Documents/UNLV/Year4/nmf3D/RcppML_nmf")

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

num_labels = 6
num_runs = 100

model <- nmf(matCSC, k = num_labels, maxit = num_runs, tol = 1e-100)

# saveRDS(model,paste0("model_",num_labels,".rds"))

# ONLY LOCALLY, NON DEV VERSION OF RcppML
colnames(model$h) <- colnames(matCSC)
rownames(model$h) <- 1:num_labels

png(paste0("h_",num_labels,".png"),width = 1300, height = 1000, res = 200)
pheatmap(t(model$h), scale = "row", border_color = NA, cluster_rows = F,
         cluster_cols = T, fontsize = 14,
         main = paste0("H matrix k=",num_labels),
         breaks = seq(-1, 1, length.out = 100))
# main = "Segmentation of H1 and naive hESCs",
# breaks = seq(-1, 1, length.out = 100)
dev.off()

# # ONLY LOCALLY, NON DEV VERSION OF RcppML
colnames(model$w) <- 1:num_labels
rownames(model$w) <- rownames(matCSC)


z <- as.data.frame(ifelse(grepl("endo",rownames(matCSC)),"endo","h1"))
rownames(z) <- rownames(matCSC)
colnames(z) <- "cell_type"

# small <- model$w[c(1:50,20000:20050),]

png(paste0("w_",num_labels,".png"),width = 1300, height = 1000, res = 200)
pheatmap(model$w, scale = "row", border_color = NA, cluster_rows = T,
         cluster_cols = F, fontsize = 14, annotation_row = z,show_rownames = F,
         main = paste0("W matrix k=",num_labels))
dev.off()




# GET LABELS:

con_mat <- con_mat %>% separate(region,c("chr","type","pos1","pos2"),sep = '[:_-]')

all_labelled <- data.frame(type = con_mat$type, chr = con_mat$chr, pos1 = con_mat$pos1, pos2 = con_mat$pos2,
                           label = colnames(model$w)[apply(model$w,1,which.max)])

h1_labelled <- subset(all_labelled, type == 'h1')[,2:5]
endo_labelled <- subset(all_labelled, type == 'endo')[,2:5]

write_tsv(h1_labelled, "h1_labelled.bed", col_names = F)
write_tsv(endo_labelled, "endo_labelled.bed", col_names = F)





