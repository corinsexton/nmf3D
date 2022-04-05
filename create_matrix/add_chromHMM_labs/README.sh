#!/bin/bash


bedtools intersect -a ../merged.bed -b ../../get_contact_labels/chromHMM_coordonly.bed -F 0.3 -wa -wb > chromHMM_merged_loops.bed
