#!/bin/bash

# Step 1: get nmf labels for each region (within RcppML_nmf/run_nmf.R)

# Step 2 Intersect regions for different cell types
bedtools intersect -a ../../RcppML_nmf/h1_labelled.bed -b ../../RcppML_nmf/endo_labelled.bed -f 0.25 -F 0.25 -wa -wb > network.bed
echo -e "H1\tEndo" > network.tsv
awk -F '\t' '{print $4"\t"$8}' network.bed >> network.tsv

# Step 3 plot_sankey.R
