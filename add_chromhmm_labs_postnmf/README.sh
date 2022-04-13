#!/bin/bash

awk '{ print gensub(/(chr[0-9XY]+)/, "\\1_endo", "g");}' ../RcppML_nmf/endo_labelled.bed  > nmf_labelled.bed
awk '{ print gensub(/(chr[0-9XY]+)/, "\\1_h1", "g");}' ../RcppML_nmf/h1_labelled.bed  >> nmf_labelled.bed

bedtools intersect -a nmf_labelled.bed -b ../get_contact_labels/chromHMM_coordonly.bed -F 0.99 -wa -wb > chromHMM_merged_loops.bed

cut -f4 chromHMM_merged_loops.bed | sort | uniq > nmf_labs.txt
cut -f8 chromHMM_merged_loops.bed | sort | uniq > chromhmm_labs.txt

./bed_to_matrix_add_chromhmm.py
