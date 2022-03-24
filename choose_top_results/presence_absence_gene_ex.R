
setwd("~/Documents/UNLV/Year4/nmf3D/choose_top_results/")

library(tidyverse)


total_int <- read_tsv("nmf3D_network_all.bed", col_names = F, col_types = 'ccccccccc')

H1_present <- total_int[total_int$X8 == '.',1:4]
endo_present <- total_int[total_int$X4 == '.',5:8]

colnames(H1_present) <- c("chr", "pos1", "pos2", "label")
colnames(endo_present) <- c("chr", "pos1", "pos2", "label")


# RAW COUNT NORMALIZATION
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



## INCLUDE GENE OVERLAPS

library(GenomicRanges)
library(AnnotationHub)

# ah <- AnnotationHub()
# gc_db <- ah[["AH75120"]]

get_overlaps <- function(df, gene_exp) {
  
  if (!exists("gc_db")) {
    
    ah <- AnnotationHub()
    # query(ah, c("Gencode", "gff", "human","GRCh38","basic"))
    # gc <- ah[["AH75120"]]
    # gc <- ah[["AH49556"]]
    
    gc_db <<- ah[["AH75120"]]
  }
  
  reg <- GRanges(df$chr, 
                    IRanges(as.numeric(df$pos1), 
                            as.numeric(df$pos2)))
  
  overlaps <- findOverlaps(reg, gc_db)
  
  genes <- extractList(gc_db$gene_name, as(overlaps, "List"))
  genes <- unstrsplit(unique(genes), ";") # Needed in case more than one gene overlaps.
  annot_regions <- cbind.data.frame(df, genes) %>%
    separate_rows(genes, sep = ";") %>% subset(genes != '')
  
  merged <- merge(annot_regions,gene_exp,by.x = "genes", by.y = "gene",all.x = T)
  
  diffs <- merged %>% select(mean_H1, mean_DEC) %>%
    mutate(diff =  mean_H1 - mean_DEC) %>% subset(!is.na(mean_H1))
  

  diffs
}
 


enh_contacting_regions_H1 <- H1_present[grepl("[12345]",H1_present$label),]
enh_contacting_regions_endo <- endo_present[grepl("[12345]",endo_present$label),]



diffs_H1 <- get_overlaps(enh_contacting_regions_H1, normalized_counts)
diffs_endo <- get_overlaps(enh_contacting_regions_endo,normalized_counts)

t.test(diffs_endo$diff, diffs_H1$diff, paired = F)



z <- cbind.data.frame(cell = c(rep('H1',nrow(diffs_H1)),rep('endo',nrow(diffs_endo))),
                      diffs = c(diffs_H1$diff,diffs_endo$diff))

ggplot(z, aes(diffs, fill = cell)) +
  geom_histogram(binwidth = 50) +
  xlim(c(-2000,2000)) + ylim(0,100)



t.test(diffs_endo$mean_DEC,diffs_endo$mean_H1, paired = T,alternative = 'g')
t.test(diffs_h1$mean_DEC,diffs_h1$mean_H1, paired = T,alternative = 'g')
t.test(normalized_counts$mean_DEC,normalized_counts$mean_H1, paired = T, alternative = 'g')




# 
# 
# z <- cbind.data.frame(cell = c(rep('H1',nrow(diffs_h1)),rep('endo',nrow(diffs_endo))),
#                       diffs = c(diffs_h1$diff_h1,diffs_endo$diff_endo))
# 
# ggplot(z, aes(diffs, fill = cell)) +
#   geom_histogram(binwidth = 50) +
#   xlim(c(-2000,2000))
# 
# 
# t.test(diffs_h1$diff_h1,diffs_endo$diff_endo, paired = F)
# 
# # z <- annot_regions[annot_regions$genes != '',]


