
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


## INCLUDE GENE EXPRESSION

library(GenomicRanges)
library(AnnotationHub)

### FIX THIS PART! -- need to use regions that are in contact here, not base region.
r <- GRanges(ss$chr, IRanges(ss$pos1,ss$pos2))

ah <- AnnotationHub()
#query(ah, c("Gencode", "gff", "human","GRCh38","basic"))
gc <- ah[["AH75120"]]

overlaps <- findOverlaps(r, gc)

genes <- extractList(gc$gene_name, as(overlaps, "List"))
genes <- unstrsplit(unique(genes), ";") # Needed in case more than one gene overlaps.
annot_regions <- cbind.data.frame(ss, genes) %>% separate_rows(genes, sep = ";")

exp_data <- read_csv("../expression/GSE75748_bulk_cell_type_ec.csv")

merged <- merge(annot_regions,exp_data,by.x = "genes", by.y = "gene")

changes <- merged %>% select(genes, nmf3D_change,chromHMM_change,H1_rep1,H1_rep2,H1_rep3,H1_rep4,DEC_rep1,DEC_rep2)

changes <- changes %>% rowwise() %>% mutate(mean_H1 = mean(c(H1_rep1,H1_rep2,H1_rep3,H1_rep4)),
                              mean_DEC =mean(c(DEC_rep1,DEC_rep2)))

# z <- annot_regions[annot_regions$genes != '',]



