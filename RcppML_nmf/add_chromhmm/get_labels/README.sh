#!/bin/bash

grep 'endo' nmf4_regions.tsv | cut -f 2,3,4 > endo_nmf4_regions_hg38.bed
grep 'h1' nmf4_regions.tsv | cut -f 2,3,4 > h1_nmf4_regions_hg38.bed
../../choose_top_results/util/liftOver h1_nmf4_regions_hg38.bed ../../choose_top_results/util/hg38ToHg19.over.chain.gz h1_nmf4_regions_hg19.bed x
../../choose_top_results/util/liftOver endo_nmf4_regions_hg38.bed ../../choose_top_results/util/hg38ToHg19.over.chain.gz endo_nmf4_regions_hg19.bed x
