#!/usr/bin/env python

classDict = dict()
counter = 0

ll = []

for line in open('class_meta.txt'):
	classDict[line.strip()] = counter
	ll.append(line.strip())
	counter += 1

outfile = open("matrix_connections.tsv",'w')
outfile.write("region\t" + '\t'.join(ll) + '\n')

bedfile = open('merged.bed','r')

#chr1    905000    906000    0.00391039668    TssBiv
#chr1    924000    925000    0.003736011843,0.005833890574,0.003736011843,0.009569902416    ReprPC,EnhBiv,TssFlnk,TssBiv
#chr1    940000    941000    0.005833890574    TssBiv
#chr1    942000    943000    0.00391039668,0.00391039668,0.00391039668    TssFlnk,ReprPC,TssBiv

class_len = len(ll)

for line in bedfile:

	outlist = ['0'] * class_len

	ll = line.strip().split()

	scores = ll[3].split(',')
	labs = ll[4].split(',')

	for i in range(len(scores)):
		outlist[classDict[labs[i]]] = str(scores[i])

	outfile.write(ll[0] + ':' + ll[1] + '-' + ll[2] + '\t' + '\t'.join(outlist) +'\n')
		


bedfile.close()
outfile.close()
