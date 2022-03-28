
setwd("~/Documents/UNLV/Year4/nmf3D/choose_top_results/")

library(tidyverse)
library(ggsignif)


total_int <- read_tsv("nmf3D_network_all.bed", col_names = F, col_types = 'ccccccccc')

H1_present <- total_int[total_int$X4 != '.',1:4]
H1_only <- total_int[total_int$X8 == '.',1:4]
endo_only <- total_int[total_int$X4 == '.',5:8]
endo_present <- total_int[total_int$X8 != '.',5:8]

colnames(H1_present) <- c("chr", "pos1", "pos2", "label")
colnames(endo_present) <- c("chr", "pos1", "pos2", "label")

colnames(H1_only) <- c("chr", "pos1", "pos2", "label")
colnames(endo_only) <- c("chr", "pos1", "pos2", "label")

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
 


enh_contacting_regions_H1 <- H1_present[grepl("[5]",H1_present$label),]
enh_contacting_regions_endo <- endo_only[grepl("[123456]",endo_only$label),]


diffs_H1 <- get_overlaps(enh_contacting_regions_H1, normalized_counts)
diffs_endo <- get_overlaps(enh_contacting_regions_endo,normalized_counts)


# h1 non labelled vs h1 label 5 expression
t.test(diffs_H1$mean_H1, diffs_endo$mean_H1, paired = F, alternative = 'g')
wilcox.test(diffs_H1$mean_H1, diffs_endo$mean_H1, paired =  F ,alternative = 'g')


H1_one <- get_overlaps(H1_present[grepl("[1]",H1_present$label),],normalized_counts)
H1_two <- get_overlaps(H1_present[grepl("[2]",H1_present$label),],normalized_counts)
H1_three <- get_overlaps(H1_present[grepl("[3]",H1_present$label),],normalized_counts)
H1_four <- get_overlaps(H1_present[grepl("[4]",H1_present$label),],normalized_counts)
H1_five <- get_overlaps(H1_present[grepl("[5]",H1_present$label),],normalized_counts)
H1_six <- get_overlaps(H1_present[grepl("[6]",H1_present$label),],normalized_counts)


z <- cbind.data.frame(cell = c(rep('label_1',nrow(H1_one)),
                               rep('label_2',nrow(H1_two)),
                               rep('label_3',nrow(H1_three)),
                               rep('label_4',nrow(H1_four)),
                               rep('label_5',nrow(H1_five)),
                               rep('label_6',nrow(H1_six)),
                               rep('no_label',nrow(diffs_endo))),
                      diffs = c(H1_one$mean_H1,
                                H1_two$mean_H1,
                                H1_three$mean_H1,
                                H1_four$mean_H1,
                                H1_five$mean_H1,
                                H1_six$mean_H1,
                                diffs_endo$mean_H1))


png("gene_exp_counts.png",width = 1300, height = 800, res = 200)
ggplot(z, aes(y = log(diffs+1), x = cell)) +
  geom_boxplot() +
  labs(x = "Regions", y = "log(normalized counts)",
       title = "Gene Expression (at base region) distribution by label") +
  geom_signif(
    comparisons = list(c("label_1","no_label"),
                        c("label_2","no_label"),
                        c("label_3","no_label"),
                        c("label_4","no_label"),
                        c("label_5","no_label"),
                        c("label_6","no_label")),
    test = 'wilcox.test',
    step_increase = c(.1),
    map_signif_level = F) + 
  theme_minimal()

dev.off()

##### GET CONTACT BOXPLOTS

contact_bed <- read_tsv("../get_contact_labels/connections.tsv",col_names = c("chr","pos1","pos2","chr_merged","pos1_merged","pos2_merged","lab_merged"),
                        col_types = "cddcccc")

sep_contact_bed <- contact_bed %>% separate_rows(chr_merged,pos1_merged,pos2_merged,sep = ',') %>% distinct()

colnames(sep_contact_bed) <- c("chr","pos1","pos2","contact_chr", "contact_pos1","contact_pos2",'labels')

endo_contacts <- sep_contact_bed[grep("_endo", sep_contact_bed$chr),]
h1_contacts <- sep_contact_bed[grep("_h1", sep_contact_bed$chr),]


h1_contacts$contact_chr <- gsub("_.*",'',h1_contacts$contact_chr)
endo_contacts$contact_chr <- gsub("_.*",'',endo_contacts$contact_chr)

h1_contacts$chr <- gsub("_.*",'',h1_contacts$chr)
endo_contacts$chr <- gsub("_.*",'',endo_contacts$chr)



#### INTERSECT nmf3D labels



nmf3D_regions <- GRanges(H1_present$chr, IRanges(as.numeric(H1_present$pos1),
                                                 as.numeric(H1_present$pos2)))
h1_regions <- GRanges(h1_contacts$chr, IRanges(as.numeric(h1_contacts$pos1),
                                                   as.numeric(h1_contacts$pos2)))



endo_only_regions <- GRanges(endo_only$chr, IRanges(as.numeric(endo_only$pos1),
                                                 as.numeric(endo_only$pos2)))
endo_regions <- GRanges(endo_contacts$chr, IRanges(as.numeric(endo_contacts$pos1),
                                               as.numeric(endo_contacts$pos2)))


x <- findOverlaps(nmf3D_regions,h1_regions)
y <- findOverlaps(endo_only_regions,endo_regions)


h1_matches <- cbind.data.frame(h1_contacts[attr(x,'to'),],H1_present[attr(x,'from'),])
h1_matches <- h1_matches %>% select(contact_chr,contact_pos1,contact_pos2, label)
colnames(h1_matches) <- c("chr", "pos1","pos2","label")

endo_matches <- cbind.data.frame(endo_contacts[attr(x,'to'),],endo_only[attr(x,'from'),])
endo_matches <- endo_matches %>% select(contact_chr,contact_pos1,contact_pos2, label)
colnames(endo_matches) <- c("chr", "pos1","pos2","label")


H1_one <- get_overlaps(h1_matches[grepl("[1]",h1_matches$label),],normalized_counts)
H1_two <- get_overlaps(h1_matches[grepl("[2]",h1_matches$label),],normalized_counts)
H1_three <- get_overlaps(h1_matches[grepl("[3]",h1_matches$label),],normalized_counts)
H1_four <- get_overlaps(h1_matches[grepl("[4]",h1_matches$label),],normalized_counts)
H1_five <- get_overlaps(h1_matches[grepl("[5]",h1_matches$label),],normalized_counts)
H1_six <- get_overlaps(h1_matches[grepl("[6]",h1_matches$label),],normalized_counts)

no_label <- get_overlaps(h1_matches[grepl("[123456]",endo_matches$label),],normalized_counts)


z <- cbind.data.frame(cell = c(rep('label_1',nrow(H1_one)),
                               rep('label_2',nrow(H1_two)),
                               rep('label_3',nrow(H1_three)),
                               rep('label_4',nrow(H1_four)),
                               rep('label_5',nrow(H1_five)),
                               rep('label_6',nrow(H1_six)),
                               rep('no_label',nrow(no_label))),
                      diffs = c(H1_one$mean_H1,
                                H1_two$mean_H1,
                                H1_three$mean_H1,
                                H1_four$mean_H1,
                                H1_five$mean_H1,
                                H1_six$mean_H1,
                                no_label$mean_H1))


png("gene_exp_counts_contacts.png",width = 1300, height = 800, res = 200)
ggplot(z, aes(y = log(diffs+1), x = cell)) +
  geom_boxplot() +
  labs(x = "Regions", y = "log(normalized counts)",
       title = "Gene Expression (at contact region) distribution by label") +
  geom_signif(
    comparisons = list(c("label_1","no_label"),
                       c("label_2","no_label"),
                       c("label_3","no_label"),
                       c("label_4","no_label"),
                       c("label_5","no_label"),
                       c("label_6","no_label")),
    test = 'wilcox.test',
    step_increase = c(.1),
    map_signif_level = F) + 
  theme_minimal()

dev.off()

table(z$cell)
