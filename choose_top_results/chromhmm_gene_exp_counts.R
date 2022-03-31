
setwd("~/Documents/UNLV/Year4/nmf3D/choose_top_results/")

library(tidyverse)
library(ggsignif)


H1_chromhmm <- read_tsv("../chromHMM_epimap_calls/H1.bed", col_names = F)
H1_chromhmm <- H1_chromhmm[,1:4]

endo_chromhmm <- read_tsv("../chromHMM_epimap_calls/endoderm.bed", col_names = F)
endo_chromhmm <- endo_chromhmm[,1:4]

colnames(H1_chromhmm) <-  c('chr','pos1','pos2','label')
colnames(endo_chromhmm) <-  c('chr','pos1','pos2','label')

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
    mutate(diff =  mean_H1 - mean_DEC) %>% subset(!is.na(mean_H1)) %>% unique()
  
  
  diffs
}


x <- unique(H1_chromhmm$label)

datalist = list()

for (i in x) {
  genes <- get_overlaps(H1_chromhmm[grepl(i,H1_chromhmm$label),],normalized_counts)
  z = data.frame(cell=c(rep(i,nrow(genes))),diffs=genes$mean_H1)
  datalist[[i]] <- z
}

all <- bind_rows(datalist)


# png("chromhmm_gene_exp_H1.png",width = 1300, height = 800, res = 200)
h1_plot <- ggplot(all, aes(y = log(diffs+1), x = cell)) +
  geom_boxplot() +
  labs(x = "Regions", y = "log(normalized counts)",
       title = "H1 Gene Expression (at base region) distribution by ChromHMM label") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# dev.off()


x <- unique(endo_chromhmm$label)

datalist = list()

for (i in x) {
  genes <- get_overlaps(endo_chromhmm[grepl(i,endo_chromhmm$label),],normalized_counts)
  z = data.frame(cell=c(rep(i,nrow(genes))),diffs=genes$mean_DEC)
  datalist[[i]] <- z
}

all <- bind_rows(datalist)

endo_plot <- ggplot(all, aes(y = log(diffs+1), x = cell)) +
  geom_boxplot() +
  labs(x = "Regions", y = "log(normalized counts)",
       title = "Endo Gene Expression (at base region) distribution by ChromHMM label") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

library(egg)


png("chromhmm_gene_exp.png",width = 1600, height = 1800, res = 200)
grid.arrange(h1_plot,endo_plot,nrow = 2,
             top = "ChromHMM gene expression")
dev.off()
 
