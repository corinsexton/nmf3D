#!/bin/bash

# Step 1: get nmf labels for each region (within RcppML_nmf/run_nmf.R)

# Step 2 Intersect regions for different cell types
bedtools intersect -a ../../RcppML_nmf/h1_labelled.bed -b ../../RcppML_nmf/endo_labelled.bed -wao  > network_h1.bed
bedtools intersect -b ../../RcppML_nmf/h1_labelled.bed -a ../../RcppML_nmf/endo_labelled.bed -wao  > network_endo.bed

cat network_h1.bed > network_all.bed

#chrX	18405000	18410000	1	.	-1	-1	.	0
#chrX	18890000	18895000	1	.	-1	-1	.	0
#chr10   785000  790000  5       chr10   785000  790000  3       5000


cat network_endo.bed | grep '\-1' - | awk '{print ".\t-1\t-1\t.\t"$1"\t"$2"\t"$3"\t"$4"\t0"};' network_endo.bed >> network_all.bed  



bedtools intersect -a ../../RcppML_nmf/h1_labelled.bed -b ../../RcppML_nmf/endo_labelled.bed -wa -wb -f .99 -F .99  > network.bed

echo -e "H1\tEndo" > network.tsv
awk -F '\t' '{print $4"\t"$8}' network.bed >> network.tsv

# Step 3 plot_sankey.R
