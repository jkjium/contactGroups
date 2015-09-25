# -*- coding: utf-8 -*-
"""
Created on Thu Sep 23 17:20:55 2015

read aligned clusters FASTA
output conserved clusters

@author: kjia
"""

# -*- coding: utf-8 -*-
import sys
import operator
from copy import deepcopy
from itertools import groupby

def fasta_iter(fasta_name):
    """
    given a fasta file. yield tuples of header, sequence
    """
    fh = open(fasta_name)
    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:

        # drop the ">"
        header = header.next()[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.next())
        yield header, seq

# return index map between group str and msa str
# key is msa idx and value is str idx
def getIdxMap(pdbseq, msaseq):
	idxMap={}
	token = -1
	for i in xrange(0,len(pdbseq)):
		for j in xrange(token+1,len(msaseq)):
			if pdbseq[i]==msaseq[j]:
				idxMap[j]=i
				token = j
				break
	if len(idxMap)!=len(pdbseq):
		print "getIdxMap()::error:length not match"
		print len(idxMap)
		print len(pdbseq)
	            
	return idxMap


def groupByCols(fa, cols):
	# example of fa: 1aba, 7.51, DGKNPRS, ++-HHPP,9 13 15 16 18 19 22
	groupArray = fa[0].split(',')
	groupStr = groupArray[2]
	groupTypeStr = groupArray[3]
	groupIdxArray = groupArray[4].split(' ')

	#print groupStr
	#print groupTypeStr
	#print groupIdxArray[0] + ' ' + groupIdxArray[1]

	idxmap = getIdxMap(groupStr, fa[1])
	#print fa[0]
	#for key in idxmap:
	#	print '%d - %d' % (key, idxmap[key])

	#for i in cols:
	#	print i

	newGroupStr = ''
	newGroupTypeStr = ''
	newGroupIdx = ''

	for i in cols:
		if i in idxmap.keys():
			str_idx = idxmap[i]
			newGroupStr = newGroupStr + groupStr[str_idx]
			newGroupTypeStr = newGroupTypeStr + groupTypeStr[str_idx]
			newGroupIdx = newGroupIdx + groupIdxArray[str_idx] + ' '

	return '%s,%s,%s,%s,%s\n' % (groupArray[0], groupArray[1], newGroupStr, newGroupTypeStr, newGroupIdx.lstrip())




def main():

	print 'enter'
	if len(sys.argv)< 4:
		print "Usage: proc_findConserved.py clusters_aligned.fa clusters_conserved.txt cluster_size"
		return 
  
	infile = sys.argv[1]
  	outfile = sys.argv[2]
  	cluster_size = int(sys.argv[3])

 	print 'infile: %s\noutfile: %s\ncluster_size: %d\n' % (infile, outfile, cluster_size)

 	print 'reading fa ...'
 	fa=fasta_iter(infile)

 	conserveDict = {}
 	aligned_fa = [] # array to hold fa records
 	# check the columns that have top number characters
 	for s in fa:
 		cluster = s[1]
 		aligned_fa.append(s)
 		for i in xrange(0,len(cluster)):
 			if cluster[i]!='-':
 				conserveDict[i] = conserveDict.get(i, 0) + 1

 	# sort the conserve Dict as decending order, return as a list of tuples
	sorted_conserve = sorted(conserveDict.items(), key=operator.itemgetter(1), reverse=True)	# 0 would be sorting the key

	cols = []
	# get top numbered columns
	for count in xrange(0,cluster_size):
		column_n_freq =  sorted_conserve[count] 
		cols.append(column_n_freq[0])
	cols.sort() # to get sorted output

	print 'writing conserved groups to %s...' % outfile
	fout = open(outfile, 'w')
	for c in aligned_fa:
	#	print c
		fout.write(groupByCols(c,cols))



if __name__=="__main__":
	main()    