#!/usr/bin/env python

classDict = dict()
counter = 0

ll = []
ll2 = []

for line in open('../class_meta.txt'):
	classDict[line.strip()] = counter
	ll.append('contact_' + line.strip())
	ll2.append(line.strip())
	counter += 1

outfile = open("matrix_connections_add_chromhmm.tsv",'w')
outfile.write("region\t" + '\t'.join(ll) + '\t' + '\t'.join(ll2) + '\n')

bedfile = open('chromHMM_merged_loops.bed','r')

# chr10_endo    510000    515000    0.01357866557,0.01357866557,0.01357866557,0.01357866557,0.01357866557,0.01357866557,0.02715733115    ReprPCWk,EnhA1,ReprPC,EnhBiv,TssBiv,TssFlnkU,EnhWk    chr10_endo    510660    511060    EnhBiv


class_len = len(ll)

for line in bedfile:

	outlist = ['0'] * class_len
	outlist2 = ['0'] * class_len

	ll = line.strip().split()

	scores = ll[3].split(',')
	labs = ll[4].split(',')

	base_lab = ll[8]
	base_lab_index = classDict[base_lab]

	for i in range(len(scores)):
		outlist[classDict[labs[i]]] = str(scores[i])
		if i == base_lab_index:
			outlist2[i] = '1'


	outfile.write(ll[5] + ':' + ll[6] + '-' + ll[7] + '\t' + '\t'.join(outlist) + '\t' + '\t'.join(outlist2) + '\n')
		


bedfile.close()
outfile.close()
