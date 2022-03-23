
setwd("~/Documents/UNLV/Year4/nmf3D/choose_top_results/")

library(tidyverse)


total_int <- read_tsv("nmf3D_network_all.bed", col_names = F, col_types = 'ccccccccc')

H1_present <- total_int[total_int$X8 == '.',1:4]
endo_present <- total_int[total_int$X4 == '.',5:8]

colnames(H1_present) <- c("chr", "pos1", "pos2", "label")
colnames(endo_present) <- c("chr", "pos1", "pos2", "label")


enh_contacting_regions_H1 <- H1_present[H1_present$label == '3',]
# enh_contacting_regions_H1 <- endo_present[endo_present$label == '5',]



## INCLUDE GENE EXPRESSION

library(GenomicRanges)
library(AnnotationHub)

ah <- AnnotationHub()
# query(ah, c("Gencode", "gff", "human","GRCh38","basic"))
gc <- ah[["AH75120"]]
# gc <- ah[["AH49556"]]

h1_reg <- GRanges(enh_contacting_regions_H1$chr, 
                  IRanges(as.numeric(enh_contacting_regions_H1$pos1), 
                          as.numeric(enh_contacting_regions_H1$pos2)))

overlaps <- findOverlaps(h1_reg, gc)

genes <- extractList(gc$gene_name, as(overlaps, "List"))
H1_genes <- unstrsplit(unique(genes), ";") # Needed in case more than one gene overlaps.
annot_regions <- cbind.data.frame(enh_contacting_regions_H1, H1_genes) %>%
  separate_rows(H1_genes, sep = ";")

# endo_reg <- GRanges(annot_regions$endo_chr, IRanges(as.numeric(annot_regions$endo_pos1),as.numeric(annot_regions$endo_pos2)))
# 
# overlaps <- findOverlaps(endo_reg, gc)
# 
# genes <- extractList(gc$gene_name, as(overlaps, "List"))
# endo_genes <- unstrsplit(unique(genes), ";") # Needed in case more than one gene overlaps.
# annot_regions <- cbind.data.frame(annot_regions, endo_genes) %>% separate_rows(endo_genes, sep = ";")

exp_data_raw <- read_csv("../expression/GSE75748_bulk_cell_type_ec.csv")
exp_data <- round(exp_data_raw[,c(2,3,4,5,9,10)])

meta <- cbind.data.frame(cell_type = c(rep("H1",4),"DEC","DEC"))
rownames(meta) <- colnames(exp_data)

library(DESeq2)

dds <- DESeqDataSetFromMatrix(countData = exp_data, colData = meta, design = ~ cell_type)
dds <- estimateSizeFactors(dds)
normalized_counts <- data.frame(counts(dds, normalized=TRUE))

normalized_counts$gene <- exp_data_raw$gene
normalized_counts <- normalized_counts %>% rowwise() %>% 
  mutate(mean_H1 = mean(c(H1_rep1,H1_rep2,H1_rep3,H1_rep4)),
         mean_DEC =mean(c(DEC_rep1,DEC_rep2)))




merged_h1 <- merge(annot_regions,normalized_counts,by.x = "H1_genes", by.y = "gene",all.x = T)
# merged_endo <- merge(annot_regions,normalized_counts,by.x = "endo_genes", by.y = "gene")

changes <- merged_h1 %>% select(mean_H1, mean_DEC)

# z <- annot_regions[annot_regions$genes != '',]

t.test(changes$mean_DEC,changes$mean_H1, paired = T,alternative = 'g')
t.test(normalized_counts$mean_DEC,normalized_counts$mean_H1, paired = T, alternative = 'g')


