#!/share/apps/R/R-4.0.2/bin/Rscript


setwd("/home/csexton/compute/NMF_differentiation/run_NMF")
library(NMF)
library(tidyverse)

con_mat <- read_tsv("matrix_connections.tsv")

x <- as.matrix(con_mat[,-1])

num_labels = 6
algorithm = 'brunet'
num_runs = 100 

res.multirun <- nmf(t(x), num_labels,method = algorithm, nrun=num_runs, .opt='P4')
 saveRDS(res.multirun,paste0("res_",num_runs,"run_",num_labels,"-",algorithm,".rds"))

#res.multirun = readRDS(paste0("res_multirun_",num_labels,".rds"))

png(paste0("coef_heatmaps/coefmap_",num_runs,"run_",num_labels,'-',algorithm,".png"),width = 2000,height = 2000, res = 200)
coefmap(res.multirun)
dev.off()

png(paste0("basis_heatmaps/asismap_",num_runs,"run_",num_labels,'-',algorithm,".png"),width = 2000,height = 2000, res = 200)
basismap(res.multirun)
dev.off()

png(paste0("consensus_heatmaps/consensusmap_",num_runs,"run_",num_labels,'-',algorithm,".png"),width = 2000,height = 2000, res = 200)
consensusmap(res.multirun)
dev.off()
