# get bedfiles


# setwd("/home/csexton/compute/NMF_differentiation/RcppML_nmf")
setwd("/Users/coripenrod/Documents/UNLV/Year4/nmf3D/RcppML_nmf/nmf4_5_chromhmmadded/")

# devtools::install_github("zdebruine/RcppSparse") # compile dev version
# devtools::install_github("zdebruine/RcppML") # compile dev version
# install.packages("RcppML")

library(Matrix)
library(RcppML)
library(tidyverse)
library(pheatmap)
library(cowplot)

con_mat <- read_tsv("../../add_chromhmm_labs_postnmf/matrix_connections_add_chromhmm.tsv")
con_mat <- con_mat %>% separate(region,c("chr","type","pos1","pos2"),sep = '[:_-]')

model <- readRDS("../with_chromhmm_model_6.rds")

scaled_h <- apply(model$h, 2, function(x) x/sum(x))

all_labelled <- data.frame(type = con_mat$type,
                           chr = con_mat$chr,
                           pos1 = con_mat$pos1,
                           pos2 = con_mat$pos2,
                           top_region = apply(scaled_h,2, max),
                           label = rownames(model$h)[apply(scaled_h,2, function(x) which.max(x))])

all_labelled <- all_labelled[all_labelled$top_region > 0.75,]

h1_labelled <- subset(all_labelled, type == 'h1')[,c(2:4,6)]
endo_labelled <- subset(all_labelled, type == 'endo')[,c(2:4,6)]

write_tsv(h1_labelled, "h1_labelled_chromhmm.bed", col_names = F)
write_tsv(endo_labelled, "endo_labelled_chromhmm.bed", col_names = F)



# H1 19/169 in HACER (39/169 all)
# endo in HACER (38/169 all)
nmf4_5_regions <- all_labelled[apply(scaled_h,2, function(x) x[["nmf4"]] > 0.3 & x[["nmf5"]] > 0.3),]
write_tsv(nmf4_5_regions,"nmf4_nmf5_regions.tsv")


nmf4 <- all_labelled[apply(scaled_h,2, function(x) x[["nmf4"]] > 0.7),]
write_tsv(nmf4,"nmf4_regions.tsv")






