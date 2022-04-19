
setwd("~/Documents/UNLV/Year4/nmf3D/intersect_v_clusters/")

library(tidyverse)

just_chromhmm <- read_tsv("../chromHMM_epimap_calls/H1.bed", col_names = c("chr","pos1","pos2", 'label','x','y','z','zz','yy')) %>% subset(grepl("EnhA",label))


raw_intersects_1 <- read_tsv("intersect_loop_HMM/raw_intersects_H1_1.bed", col_names = c("chr_loop1","pos1_loop1","pos2_loop1",
                                                                                     "chr_loop2","pos1_loop2","pos2_loop2", "score",
                                                                                     "chr","pos1","pos2",'label','x','y','z','zz','yy')) %>% subset(grepl("EnhA",label))


raw_intersects_2 <- read_tsv("intersect_loop_HMM/raw_intersects_H1_2.bed", col_names = c("chr_loop1","pos1_loop1","pos2_loop1", "score",
                                                                                       "chr","pos1","pos2",'label','x','y','z','zz','yy')) %>% subset(grepl("EnhA",label))

raw_intersects <- rbind.data.frame(raw_intersects_1[,c(1,2,3,7,8,9,10,11,12,13,14,15,16)],raw_intersects_2)


nmf_intersects <- read_tsv("../choose_top_results/total_network_H1.bed") %>% subset(grepl("EnhA",chromHMM_H1) & grepl("nmf[245]",nmf3D_H1))

clusters <- read_tsv("../RcppML_nmf/h1_labelled.bed",col_names = c('chr','pos1','pos2','label')) %>% subset(grepl("nmf[245]",label))
clusters_with_chromhmm <- read_tsv("../RcppML_nmf/add_chromhmm/h1_labelled_chromhmm.bed",col_names = c('chr','pos1','pos2','label')) %>% subset(grepl("nmf[245]",label))



# ANNOTATION OF ENHANCERS FROM HACER DATABASE
enhancers <- read_tsv("../annotation/GROseq-H1_ESCs/GROseq-H1-HACER.tsv")

# ANNOTATION OF ENHANCERS FROM HACER DATABASE
enhancers <- read_tsv("../annotation/H1_enhancerAtlas.bed",col_names = c("chr","start","end",'score'))

library(GenomicRanges)
enh_regions <- GRanges(enhancers$chr,IRanges(enhancers$start,enhancers$end))

#### FIND PRESENCE OF ENHANCERS IN EACH DATASET

find_overlaps_enh <- function(bed, enhancer_regions) {
  
  query_regions <- GRanges(bed$chr,IRanges(bed$pos1,bed$pos2))
  
  overlaps <- findOverlaps(query_regions, enhancer_regions)
  
  cat(paste0("True Positive: ",round(length(overlaps) / nrow(bed),3),'\n'))
  cat(paste0("False Positive: ", round((nrow(bed) -length(overlaps)) / nrow(bed),3),'\n'))
  cat(paste0("False Negative: ", round((length(enhancer_regions) -length(overlaps)) /length(enhancer_regions),3)))
  
  
}

find_overlaps_enh(just_chromhmm,enh_regions)
find_overlaps_enh(raw_intersects,enh_regions)
find_overlaps_enh(nmf_intersects,enh_regions)
find_overlaps_enh(clusters,enh_regions)
find_overlaps_enh(clusters_with_chromhmm,enh_regions)





genes <- extractList(gc$gene_name, as(overlaps, "List"))
H1_genes <- unstrsplit(unique(genes), ";") # Needed in case more than one gene overlaps.
annot_regions <- cbind.data.frame(total_matches, H1_genes) %>% separate_rows(H1_genes, sep = ";")

