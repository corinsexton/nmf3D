#!/bin/bash

mkdir labelled
rm class_meta.txt

for bedfile in ../get_contact_labels/state_bedfiles/scored/*bed; do
	
	id=${bedfile%_scored.bed}
	id=${id##*/}
	echo $id >> class_meta.txt

	awk -v id=$id '{print $1"\t"$2"\t"$3"\t"$4"\t"id}' $bedfile > labelled/${id}.bed


done

cat labelled/* > labelled_all.bed

bedtools sort -i labelled_all.bed | bedtools merge -i - -d -1000 -o collapse -c 4,5 > merged.bed


./bed_to_matrix.py
