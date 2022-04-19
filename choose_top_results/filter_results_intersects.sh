#!/bin/bash

# Step 1: get nmf labels for each region (within RcppML_nmf/run_nmf.R)

# Step 2 Intersect regions for different cell types

 ## ALL REGIONS
 bedtools intersect -a ../RcppML_nmf/h1_labelled.bed -b ../RcppML_nmf/endo_labelled.bed -wao  > nmf3D_network_h1.bed
 bedtools intersect -b ../RcppML_nmf/h1_labelled.bed -a ../RcppML_nmf/endo_labelled.bed -wao  > nmf3D_network_endo.bed
# 
cat nmf3D_network_h1.bed > nmf3D_network_all.bed
# 
# #chrX	18405000	18410000	1	.	-1	-1	.	0
# #chrX	18890000	18895000	1	.	-1	-1	.	0
# #chr10   785000  790000  5       chr10   785000  790000  3       5000
# 
cat nmf3D_network_endo.bed | grep '\-1' - | awk '{print ".\t-1\t-1\t.\t"$1"\t"$2"\t"$3"\t"$4"\t0"};' - >> nmf3D_network_all.bed  



## JUST INTERSECTING REGIONS

# nmf3D regions

#echo -e "chr\tpos1\tpos2\tH1_nmf3D_lab\tEndo_nmf3D_lab" > nmf3D_network.tsv
bedtools intersect -a ../RcppML_nmf/h1_labelled.bed -b ../RcppML_nmf/endo_labelled.bed -wa -wb -f .99 -F .99  | \
	awk -F '\t' '{print $1"\t"$2"\t"$3"\t"$4"\t"$8}' - > nmf3D_network.bed


# chromHMM regions
#echo -e "chr\tpos1\tpos2\tH1_chromHMM_lab\tEndo_chromHMM_lab" > chromHMM_network.tsv
bedtools intersect -a ../chromHMM_epimap_calls/H1.bed -b ../chromHMM_epimap_calls/endoderm.bed -wa -wb -f .99 -F .99 | \
	awk -F '\t' '{print $1"\t"$2"\t"$3"\t"$4"\t"$13}' - | sed 's/\//_/g' > chromHMM_network.bed

# contacts
bedtools intersect -a ../chromHMM_epimap_calls/H1.bed -b ../chromHMM_epimap_calls/endoderm.bed -wa -wb -f .99 -F .99 | \
	awk -F '\t' '{print $1"\t"$2"\t"$3"\t"$4"\t"$13}' - | sed 's/\//_/g' > chromHMM_network.bed


### TOTAL INTERSECT
# a windows = 5kb
# b windows = 200 - 20kb, generally smaller

echo -e "chr\tpos1\tpos2\tnmf3D_H1\tnmf3D_endo\tchr_1\tpos1_1\tpos2_1\tchromHMM_H1\tchromHMM_endo" > total_network.bed
bedtools intersect -a  nmf3D_network.bed -b chromHMM_network.bed -wa -wb -f .05 -F .99 >> total_network.bed
	#awk -F '\t' '{print $1"\t"$2"\t"$3"\t"$4"\t"$13}' - | sed 's/\//_/g' >> chromHMM_network.tsv


#chr10	785000	790000	nmf4	.	-1	-1	.	0	chr10	785460	786060	ReprPC	0	.	785460	786060	128,128,128

echo -e "chr\tpos1\tpos2\tnmf3D_H1\tchr_chromhmm\tpos1_chromhmm\tpos2_chromhmm\tchromHMM_H1" > total_network_H1.bed
cut -f 1,2,3,4 nmf3D_network_h1.bed | bedtools intersect -a - -b ../chromHMM_epimap_calls/H1.bed -wa -wb -f .05 -F .99 | cut -f 1,2,3,4,5,6,7,8 >> total_network_H1.bed
	#awk -F '\t' '{print $1"\t"$2"\t"$3"\t"$4"\t"$13}' - | sed 's/\//_/g' >> chromHMM_network.tsv


