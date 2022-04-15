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
library(cowplot)

con_mat <- read_tsv("../add_chromhmm_labs_postnmf/matrix_connections_add_chromhmm.tsv")

x <- as.matrix(con_mat[,-1])
rownames(x) <- con_mat$region
matCSC <- as(t(x), "dgCMatrix")

# # Selecting Rank
# 
# cv_predict <- crossValidate(matCSC, k = 1:20, method = "predict", reps = 3, seed = 123)
# cv_robust <- crossValidate(matCSC, k = 1:20, method = "robust", reps = 3, seed = 123)
# cv_impute <- crossValidate(matCSC, k = 1:20, method = "impute", reps = 3, seed = 123)
# 
# png("rank_selection.png", width = 1500, height = 500, res = 200)
# plot_grid(
#   plot(cv_predict) + ggtitle("method = 'predict'") + theme(legend.position = "none"),
#   plot(cv_robust) + ggtitle("method = 'robust'") + theme(legend.position = "none"),
#   plot(cv_impute) + ggtitle("method = 'impute'") + scale_y_continuous(trans = "log10") + theme(legend.position = "none"),
#   get_legend(plot(cv_predict)), rel_widths = c(1, 1, 1, 0.4), nrow = 1, labels = "auto")
# dev.off()



num_labels = 6
num_runs = 1000

model <- nmf(matCSC, k = num_labels, maxit = num_runs, tol = 1e-20)

# saveRDS(model,paste0("../with_chromhmm_model_",num_labels,".rds"))

# small <- model$h[,c(1:3,10:30)]
z <- as.data.frame(ifelse(grepl("endo",colnames(model$h)),"endo","h1"))
rownames(z) <- colnames(model$h)
colnames(z) <- "cell_type"


# png(paste0("h_",num_labels,".png"),width = 1300, height = 1000, res = 200)
h <- pheatmap(model$h, scale = "column", border_color = NA, cluster_rows = F,
              cluster_cols = T, fontsize = 14, show_colnames = F,
              main = "H matrix", treeheight_col = 0,
              annotation_col = z,
              breaks = seq(-1, 1, length.out = 100))
# main = "Segmentation of H1 and naive hESCs",
# breaks = seq(-1, 1, length.out = 100)
# dev.off()

zz <- model$w
rownames(zz) <- c("EnhA1","EnhA2","EnhBiv","EnhG1","EnhG2","EnhWk","Het","Quies","ReprPC","ReprPCWk",
                       "TssA","TssBiv","TssFlnk","TssFlnkD","TssFlnkU","Tx","TxWk","ZNF-Rpts",
                       "1_EnhWk/Quies","2_Tx","3_Biv","4_TSS","5_Enh","6_Repr")

order_rows <- c("EnhA1","EnhA2","EnhG1","EnhG2","EnhWk","TssA","TssFlnk","TssFlnkD","TssFlnkU",
                "Tx","TxWk","EnhBiv","TssBiv","ReprPC","ReprPCWk","Het","Quies", "ZNF-Rpts",
                "1_EnhWk/Quies","2_Tx","3_Biv","4_TSS","5_Enh","6_Repr")

# png(paste0("w_",num_labels,".png"),width = 1300, height = 1000, res = 200)
w <- pheatmap(zz[order_rows,], scale = "row", border_color = NA, cluster_rows = F,
              cluster_cols = F, fontsize = 14, show_rownames = T,
              breaks = seq(-1, 1, length.out = 100),
              main = "W matrix")

# main = paste0("W matrix k=",num_labels))
# dev.off()


png(paste0("nmf3D_results_with_chromhmm",num_labels,".png"),width = 2600, height = 1000, res = 200)
plot_grid(w$gtable, h$gtable,ncol = 2, labels = "AUTO", rel_widths = c(1, 1.5))
dev.off()



# GET LABELS:

con_mat <- con_mat %>% separate(region,c("chr","type","pos1","pos2"),sep = '[:_-]')

scaled_h <- apply(model$h, 2, function(x) x/sum(x))

all_labelled <- data.frame(type = con_mat$type,
                           chr = con_mat$chr,
                           pos1 = con_mat$pos1,
                           pos2 = con_mat$pos2,
                           top_region = apply(scaled_h,2, max),
                           label = rownames(model$h)[apply(scaled_h,2, function(x) which.max(x))])

all_labelled <- all_labelled[all_labelled$top_region > 0.75,]


# H1 19/169 in HACER (39/169 all)
# endo in HACER (38/169 all)
nmf4_5_regions <- all_labelled[apply(scaled_h,2, function(x) x[["nmf4"]] > 0.3 & x[["nmf5"]] > 0.3),]
write_tsv(nmf4_5_regions,"nmf4_nmf5_regions.tsv")

h1_labelled <- subset(all_labelled, type == 'h1')[,c(2:4,6)]
endo_labelled <- subset(all_labelled, type == 'endo')[,c(2:4,6)]

write_tsv(h1_labelled, "h1_labelled_chromhmm.bed", col_names = F)
write_tsv(endo_labelled, "endo_labelled_chromhmm.bed", col_names = F)

