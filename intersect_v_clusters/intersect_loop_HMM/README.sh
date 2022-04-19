#!/bin/bash


bedtools intersect -a ../../score_loops/h1.loops.scored -b ../../chromHMM_epimap_calls/H1.bed -wa -wb -f 0.05 -F 0.99 > raw_intersects_H1_1.bed
cut -f4,5,6,7 ../../score_loops/h1.loops.scored | bedtools intersect -a -  -b ../../chromHMM_epimap_calls/H1.bed -wa -wb -f 0.05 -F 0.99 > raw_intersects_H1_2.bed
