#!/home/csexton/.conda/envs/3dgenome/bin/python

## https://cooltools.readthedocs.io/en/latest/notebooks/viz.html#Inspecting-C-data

import numpy as np
import pandas as pd
import cooltools
import cooler

# get resolutions
#cooler.fileops.list_coolers('./4DNFI9GMP2J8.mcool')

# ls ../raw_data/
# 4DNFI9E222YJ.mcool  4DNFIJLK5WML.mcool

#clr = cooler.Cooler('../raw_data/4DNFI9E222YJ.mcool::resolutions/5000')
#
#loops_file = open("../call_loops/endoderm.loops", 'r')
#outfile = open("endoderm.loops.scored",'w')

clr = cooler.Cooler('../raw_data/4DNFIJLK5WML.mcool::resolutions/5000')

loops_file = open("../call_loops/h1.loops", 'r')
outfile = open("h1.loops.scored",'w')

for line in loops_file:
	ll = line.strip().split()

	region1 = (ll[0],ll[1],ll[2])
	region2 = (ll[3],ll[4],ll[5])
	
	normalized_score = clr.matrix(balance=True).fetch(region1,region2)[0,0]

	outfile.write('\t'.join(ll[0:6]) + '\t' + str(normalized_score) + '\n')
	
	# returns score : [[0.01077166]]

outfile.close()

