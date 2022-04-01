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

con_mat <- read_tsv("../create_matrix/matrix_connections.tsv")

x <- as.matrix(con_mat[,-1])
rownames(x) <- con_mat$region
matCSC <- as(x, "dgCMatrix")

# Selecting Rank

cv_predict <- crossValidate(matCSC, k = 1:20, method = "predict", reps = 3, seed = 123)
cv_robust <- crossValidate(matCSC, k = 1:20, method = "robust", reps = 3, seed = 123)
cv_impute <- crossValidate(matCSC, k = 1:20, method = "impute", reps = 3, seed = 123)

png("rank_selection.png", width = 1500, height = 500, res = 200)
plot_grid(
  plot(cv_predict) + ggtitle("method = 'predict'") + theme(legend.position = "none"),
  plot(cv_robust) + ggtitle("method = 'robust'") + theme(legend.position = "none"),
  plot(cv_impute) + ggtitle("method = 'impute'") + scale_y_continuous(trans = "log10") + theme(legend.position = "none"),
  get_legend(plot(cv_predict)), rel_widths = c(1, 1, 1, 0.4), nrow = 1, labels = "auto")
dev.off()



num_labels = 6
num_runs = 100

model <- nmf(matCSC, k = num_labels, maxit = num_runs, tol = 1e-20)

# saveRDS(model,paste0("model_",num_labels,".rds"))

png(paste0("h_",num_labels,".png"),width = 1300, height = 1000, res = 200)
pheatmap(t(model$h), scale = "row", border_color = NA, cluster_rows = F,
         cluster_cols = T, fontsize = 14,
         main = paste0("H matrix k=",num_labels),
         breaks = seq(-2, 2, length.out = 100))
# main = "Segmentation of H1 and naive hESCs",
# breaks = seq(-1, 1, length.out = 100)
dev.off()




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

all_labelled <- data.frame(type = con_mat$type,
                           chr = con_mat$chr,
                           pos1 = con_mat$pos1,
                           pos2 = con_mat$pos2,
                           label = colnames(model$w)[apply(model$w,1, function(x) which(x>0.8))])

h1_labelled <- subset(all_labelled, type == 'h1')[,2:5]
endo_labelled <- subset(all_labelled, type == 'endo')[,2:5]

write_tsv(h1_labelled, "h1_labelled.bed", col_names = F)
write_tsv(endo_labelled, "endo_labelled.bed", col_names = F)








