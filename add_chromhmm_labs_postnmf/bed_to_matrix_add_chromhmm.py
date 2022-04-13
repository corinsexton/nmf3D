#!/usr/bin/env python

classDict = dict()
classDict2 = dict()
counter = 0

ll = []
ll2 = []

for line in open('nmf_labs.txt'):
	classDict2[line.strip()] = counter
	ll2.append(line.strip())
	counter += 1

counter = 0
for line in open('chromhmm_labs.txt'):
	classDict[line.strip()] = counter
	ll.append(line.strip())
	counter += 1

outfile = open("matrix_connections_add_chromhmm.tsv",'w')
outfile.write("region\t" + '\t'.join(ll) + '\t' + '\t'.join(ll2) + '\n')

bedfile = open('chromHMM_merged_loops.bed','r')

#chr10_endo	510000	515000	2	chr10_endo	510660	511060	EnhBiv
#chr10_endo	510000	515000	2	chr10_endo	511060	511460	EnhA2
#chr10_endo	510000	515000	2	chr10_endo	511460	511660	EnhWk
#chr10_endo	510000	515000	2	chr10_endo	511660	517060	TxWk
#chr10_endo	785000	790000	5	chr10_endo	784460	785260	TssBiv

class_len = len(ll)
class_len2 = len(ll2)

for line in bedfile:

	outlist = ['0'] * class_len
	outlist2 = ['0'] * class_len2

	ll = line.strip().split()

	nmf_lab = ll[3].strip('nmf')
	base_lab = ll[7]

	base_lab_index = classDict[base_lab]

	outlist[base_lab_index] = '1'
	outlist2[int(nmf_lab) - 1] = '1'

	outfile.write(ll[4] + ':' + ll[5] + '-' + ll[6] + '\t' + '\t'.join(outlist) + '\t' + '\t'.join(outlist2) + '\n')
		


bedfile.close()
outfile.close()
