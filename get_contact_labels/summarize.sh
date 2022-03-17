#!/bin/bash


## ADD LABELS TO ENDODERM/H1 CELLS (chr1_H1) and chromHMM

cut -f1-7 ../score_loops/endoderm.loops.scored | awk '{ print gensub(/(chr[0-9XY]+)/, "\\1_endo", "g");}' > labelled_endoderm.loops
cut -f1-7 ../score_loops/h1.loops.scored | awk '{ print gensub(/(chr[0-9XY]+)/, "\\1_h1", "g");}' > labelled_h1.loops

cat labelled_endoderm.loops labelled_h1.loops > all_loops.bedpe

awk '{ print gensub(/ZNF\/Rpts/, "ZNF-Rpts", "g");}' ../chromHMM_epimap_calls/H1.bed  > x
awk '{ print gensub(/ZNF\/Rpts/, "ZNF-Rpts", "g");}' ../chromHMM_epimap_calls/endoderm.bed  > y

awk '{ print gensub(/(chr[0-9XY]+)/, "\\1_h1", "g");}' x  > chromHMM_h1.bed
awk '{ print gensub(/(chr[0-9XY]+)/, "\\1_endo", "g");}' y > chromHMM_endoderm.bed

cat chromHMM_endoderm.bed chromHMM_h1.bed > chromHMM_total.bed


chromHMM=chromHMM_total.bed
loops=all_loops.bedpe


cut -f1-4 $chromHMM > chromHMM_coordonly.bed
cut -f4 chromHMM_coordonly.bed | sort | uniq > class_meta.txt

##### LOOPS FILE, 7th COLUMN NEEDS TO BE SCORE


### STEP 1: Intersect loops (micro-c) with chromHMM labels
bedtools intersect -a $loops -b chromHMM_coordonly.bed -wa -wb -F .3 > intersect_loops_states2.bed
awk '{print $4"\t"$5"\t"$6"\t"$1"\t"$2"\t"$3"\t"$7}' $loops > loops_coord_only1.bedpe 
bedtools intersect -a loops_coord_only1.bedpe -b chromHMM_coordonly.bed -wa -wb -F .3 > intersect_loops_states1.bed

#bedtools sort -i intersect_loops_states1.bed | bedtools merge -i - -o collapse,count,count_distinct -c 11 -d -5000 > loops1_states.bed
#bedtools sort -i intersect_loops_states2.bed | bedtools merge -i - -o collapse,count,count_distinct -c 11 -d -5000 > loops2_states.bed

bedtools sort -i intersect_loops_states1.bed | bedtools merge -i - -c 4,5,6,11 -o collapse > connections.tsv
bedtools sort -i intersect_loops_states2.bed | bedtools merge -i - -c 4,5,6,11 -o collapse >> connections.tsv

#grep '_endo' connections.tsv > connections_endo.tsv
#grep '_h1' connections.tsv > connections_h1.tsv


### STEP 2: Create 2 bedfiles (contact 1 -> 2 and contact 2 -> 1)
cut -f 4,5,6,7,11 intersect_loops_states2.bed > contact2_states.bed
cut -f 4,5,6,7,11 intersect_loops_states1.bed > contact1_states.bed



#### STEP 3: SEPARATE AND SCORE BY STATE

mkdir state_bedfiles
mkdir state_bedfiles/scored

#cut -f4 segway_run8_coordonly.bed | sort | uniq > classes_meta.txt

while read -r state; do

	grep -P "\t$state$" contact2_states.bed > state_bedfiles/${state}.bed
	grep -P "\t$state$" contact1_states.bed >> state_bedfiles/${state}.bed

	bedtools sort -i state_bedfiles/${state}.bed > z
	bedtools merge -o sum -c 4 -i z -d -5000 > state_bedfiles/scored/${state}_scored.bed

	#awk '{print $1"\t"$2"\t"$3"\t1"}' x > state_bedfiles_meta/scored/${state}_scored.bed
	#bedtools sort -i state_bedfiles_meta/scored/${state}_scored.bed > z
	#mv z state_bedfiles_meta/scored/${state}_scored.bed

	#bedtools complement -i state_bedfiles_meta/scored/${state}_scored.bed -g hg38.chrom.sizes > z
	#awk '{print $1"\t"$2"\t"$3"\t0"}' z >> state_bedfiles_meta/scored/${state}_scored.bed

	#bedtools sort -i state_bedfiles_meta/scored/${state}_scored.bed > z
	#mv z state_bedfiles_meta/scored/${state}_scored.bed


done < class_meta.txt



#### STEP 4: CREATE genome data archive for segway
#module load conda
#conda activate segway
##
##pip install genomedata
#
#
##EnhA1_scored.bed  EnhBiv_scored.bed  EnhG2_scored.bed  Het_scored.bed    ReprPC_scored.bed    TssA_scored.bed    TssFlnkD_scored.bed  TssFlnkU_scored.bed  TxWk_scored.bed
##EnhA2_scored.bed  EnhG1_scored.bed   EnhWk_scored.bed  Quies_scored.bed  ReprPCWk_scored.bed  TssBiv_scored.bed  TssFlnk_scored.bed   Tx_scored.bed        ZNF-Rpts_scored.bed
#
#
#### TOOK OVER 24 HRS (around 28 hrs on login node)
#
#genomedata-load --sizes -s hg38.chrom.sizes \
#	-t EnhA1=state_bedfiles/scored/EnhA1_scored.bed \
#	-t EnhA2=state_bedfiles/scored/EnhA2_scored.bed \
#	-t EnhBiv=state_bedfiles/scored/EnhBiv_scored.bed \
#	-t EnhG1=state_bedfiles/scored/EnhG1_scored.bed \
#	-t EnhG2=state_bedfiles/scored/EnhG2_scored.bed \
#	-t EnhWk=state_bedfiles/scored/EnhWk_scored.bed \
#	-t Het=state_bedfiles/scored/Het_scored.bed \
#	-t Quies=state_bedfiles/scored/Quies_scored.bed \
#	-t ReprPC=state_bedfiles/scored/ReprPC_scored.bed \
#	-t ReprPCWk=state_bedfiles/scored/ReprPCWk_scored.bed \
#	-t TssA=state_bedfiles/scored/TssA_scored.bed \
#	-t TssBiv=state_bedfiles/scored/TssBiv_scored.bed \
#	-t TssFlnkD=state_bedfiles/scored/TssFlnkD_scored.bed \
#	-t TssFlnk=state_bedfiles/scored/TssFlnk_scored.bed \
#	-t TssFlnkU=state_bedfiles/scored/TssFlnkU_scored.bed \
#	-t Tx=state_bedfiles/scored/Tx_scored.bed \
#	-t TxWk=state_bedfiles/scored/TxWk_scored.bed \
#	-t ZNFRpts=state_bedfiles/scored/ZNF-Rpts_scored.bed \
#	chromhmm_micro-c_1kb.data


