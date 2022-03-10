
setwd("~/Documents/UNLV/Year4/nmf3D/choose_top_results/")

library(tidyverse)


total_int <- read_tsv("total_network.bed")

table(total_int$nmf3D_H1,total_int$nmf3D_endo)
table(total_int$chromHMM_H1,total_int$chromHMM_endo)

z <- table(total_int$nmf3D_H1,total_int$nmf3D_endo,total_int$chromHMM_H1,total_int$chromHMM_endo)

total_int <- total_int %>% mutate(nmf3D_change = paste0(nmf3D_H1,'_',nmf3D_endo)) %>% mutate(chromHMM_change = paste0(chromHMM_H1,'_',chromHMM_endo))

interesting_cells <- (total_int$nmf3D_H1 == 1 | total_int$nmf3D_endo == 1) & !(total_int$nmf3D_endo == 1 & total_int$nmf3D_H1 == 1)
ss <- total_int[interesting_cells,]
z <- table(ss$nmf3D_change, ss$chromHMM_change)
zz <- data.frame(z)

nonzero <- subset(zz, Freq != 0)


write_tsv(subset(ss, nmf3D_H1 == 1 & chromHMM_endo == 'EnhA1'),"candidate_regions.tsv")

write_tsv(ss[,c(6,7,8)], "all_candidate_regions.tsv", col_names =F)
write_tsv(ss[,c(6,7,8,11,12)], "full_candidate_regions.tsv", col_names =F)
