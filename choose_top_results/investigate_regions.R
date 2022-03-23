
setwd("~/Documents/UNLV/Year4/nmf3D/choose_top_results/")

library(tidyverse)


total_int <- read_tsv("total_network.bed")

table(total_int$nmf3D_H1,total_int$nmf3D_endo)
table(total_int$chromHMM_H1,total_int$chromHMM_endo)

z <- table(total_int$nmf3D_H1,total_int$nmf3D_endo,total_int$chromHMM_H1,total_int$chromHMM_endo)

total_int <- total_int %>% mutate(nmf3D_change = paste0(nmf3D_H1,'_',nmf3D_endo)) %>% mutate(chromHMM_change = paste0(chromHMM_H1,'_',chromHMM_endo))

interesting_cells <- (total_int$nmf3D_H1 == 1 | total_int$nmf3D_endo == 1) & 
  !(total_int$nmf3D_endo == 1 & total_int$nmf3D_H1 == 1)

interesting_cells <- (grepl("Tss",total_int$chromHMM_H1,"Tss") | grepl("Tss",total_int$chromHMM_endo)) &
  !(grepl("Tss",total_int$chromHMM_H1) & grepl("Tss",total_int$chromHMM_endo))

interesting_cells <- (total_int$nmf3D_H1 == 5 | total_int$nmf3D_endo == 5)


ss <- total_int[interesting_cells,]
z <- table(ss$nmf3D_change, ss$chromHMM_change)
zz <- data.frame(z)

nonzero <- subset(zz, Freq != 0)


write_tsv(subset(ss, nmf3D_H1 == 1 & chromHMM_endo == 'EnhA1'),"candidate_regions.tsv")

write_tsv(ss[,c(6,7,8)], "all_candidate_regions.tsv", col_names =F)
write_tsv(ss[,c(6,7,8,11,12)], "full_candidate_regions.tsv", col_names =F)

# MERGE CONTACT DATA:

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

# FIND OVERLAPS IN SIGNIFICANT nmf3D LABELLED REGIONS
library(GenomicRanges)

nmf3D_regions <- GRanges(ss$chr, IRanges(ss$pos1,ss$pos2))
endo_regions <- GRanges(endo_contacts$chr, IRanges(as.numeric(endo_contacts$pos1),
                                                   as.numeric(endo_contacts$pos2)))


x <- findOverlaps(nmf3D_regions,endo_regions)

endo_matches <- cbind.data.frame(endo_contacts[attr(x,'to'),],ss[attr(x,'from'),])

h1_regions <- GRanges(h1_contacts$chr, IRanges(as.numeric(h1_contacts$pos1),as.numeric(h1_contacts$pos2)))
endo_matched_regions <- GRanges(endo_matches[,1], IRanges(as.numeric(endo_matches[,2]),as.numeric(endo_matches[,3])))

x <- findOverlaps(endo_matched_regions,h1_regions)

total_matches <- cbind.data.frame(h1_contacts[attr(x,'to'),],endo_matches[attr(x,'from'),])

colnames(total_matches) <- c("H1_chr", "H1_pos1","H1_pos2","H1_contact_chr",
                             "H1_contact_pos1", "H1_contact_pos2","H1_contact_labels",
                             "endo_chr", "endo_pos1","endo_pos2","endo_contact_chr",
                             "endo_contact_pos1", "endo_contact_pos2","endo_contact_labels",
                             "nmf3D_chr", "nmf3D_pos1", "mnf3D_pos2","nmf3D_H1_lab", "nmf3D_endo_lab", 
                             "chromHMM_chr", "chromHMM_pos1", "chromHMM_pos2",
                             "chromHMM_H1_lab", "chromHMM_endo_lab",
                             "nmf3D_change", "chromHMM_change")

## INCLUDE GENE EXPRESSION

library(GenomicRanges)
library(AnnotationHub)

# h1_reg <- GRanges(total_matches$H1_contact_chr, IRanges(as.numeric(total_matches$H1_contact_pos1),
#                                                         as.numeric(total_matches$H1_contact_pos2)))
h1_reg <- GRanges(total_matches$H1_chr, IRanges(as.numeric(total_matches$H1_pos1),
                                                        as.numeric(total_matches$H1_pos2)))

ah <- AnnotationHub()
# query(ah, c("Gencode", "gff", "human","GRCh38","basic"))
gc <- ah[["AH75120"]]
# gc <- ah[["AH49556"]]

overlaps <- findOverlaps(h1_reg, gc)

genes <- extractList(gc$gene_name, as(overlaps, "List"))
H1_genes <- unstrsplit(unique(genes), ";") # Needed in case more than one gene overlaps.
annot_regions <- cbind.data.frame(total_matches, H1_genes) %>% separate_rows(H1_genes, sep = ";")

# endo_reg <- GRanges(annot_regions$endo_contact_chr, IRanges(as.numeric(annot_regions$endo_contact_pos1),as.numeric(annot_regions$endo_contact_pos2)))

endo_reg <- GRanges(annot_regions$endo_chr, IRanges(as.numeric(annot_regions$endo_pos1),as.numeric(annot_regions$endo_pos2)))

overlaps <- findOverlaps(endo_reg, gc)

genes <- extractList(gc$gene_name, as(overlaps, "List"))
endo_genes <- unstrsplit(unique(genes), ";") # Needed in case more than one gene overlaps.
annot_regions <- cbind.data.frame(annot_regions, endo_genes) %>% separate_rows(endo_genes, sep = ";")

exp_data <- read_csv("../expression/GSE75748_bulk_cell_type_ec.csv")

merged_h1 <- merge(annot_regions,exp_data,by.x = "H1_genes", by.y = "gene",all.x = T)
merged_endo <- merge(annot_regions,exp_data,by.x = "endo_genes", by.y = "gene")

changes <- merged_h1 %>% select(H1_genes, endo_genes, nmf3D_change, chromHMM_change, H1_rep1,H1_rep2,H1_rep3,H1_rep4,DEC_rep1,DEC_rep2)

changes <- changes %>% rowwise() %>% mutate(mean_H1 = mean(c(H1_rep1,H1_rep2,H1_rep3,H1_rep4)),
                              mean_DEC =mean(c(DEC_rep1,DEC_rep2)))

# z <- annot_regions[annot_regions$genes != '',]



