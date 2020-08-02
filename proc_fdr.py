import numpy as np
from itertools import groupby

import commp as cp
from utils_pfammsa import pfammsa

def topfdrdca(args):
	if len(args)!=2:
		cp._err('Usage: python proc_fdr.py topfdrdca infile outfile')
	infile = args[0]
	outfile = args[1]
	lines = [line.split(' ') for line in cp.loadlines(infile)]
	fout = open(outfile, 'w')
	#t = [list(x[1]) for x in groupby(lines, lambda line: (line[0],line[1],line[2]))]
	# x[1]
	# [['375', '488', '521', '375', '488', '-0.354114', '0.424778'], ['375', '488', '521', '375', '521', '-0.306549', '0.796030'], ['375', '488', '521', '488', '521', '0.051714', '0.707303']]
	for x in groupby(lines, lambda line: ('%s %s %s' % (line[0],line[1],line[2]))):
		p = list(x[1])
		p1 = p[0]
		p2 = p[1]
		p3 = p[2]
		outstr = '%s %s %s %s %s %s %s %s %s %s %.8f' % (x[0], p1[3],p1[4],p1[5],p2[3],p2[4],p2[5],p3[3],p3[4],p3[5], float(p1[5])+float(p2[5])+float(p3[5]))
		fout.write('%s\n' % outstr)
	fout.close()



# convert each column of a MSA (n sequences) into a n by {1 | 2 | ... | 20} matrix
# each column matrix is saved separately using the given prefix with its column number (0-based)
# before implementing wen's filter method 
# use: 
# translate score file back to a compact msa
# $ python utils_pfammsa.py score2msa PF00043_p90.txt.score aa PF00043 PF00043.r.fa
def msa2vectors(args):
	'''
	def _expand_alphabet(scheme, ch):
		return [scheme[k] if k==ch else 1 for k in scheme]
	'''

	assert len(args) == 3, 'Usage: python proc_fdr.py msa2vectors score_reversed_fa_file .rcol output_prefix'
	infile = args[0]
	cols = ['%d' % col for col in np.loadtxt(args[1],delimiter=',')]
	outprefix = args[2]

	pfm = pfammsa(infile) # default opt='ambaa' converts all abnormal alphabets into gap '.'
	for c in xrange(0, pfm.msalen):
		print c
		clist = pfm.msacol(c)
		# convert n x column alphabets into n raw matrix
		cmatrix = [[0 if aa!=a else 1 for aa in cp.aas02] for a in clist]
		#print clist

		# convert to np array for slicing
		npcmatrix = np.array(cmatrix)
		#print npcmatrix
		# remove all-zero columns from the matrix
		#print sum(npcmatrix)
		csum = sum(npcmatrix)
		nonzero_idx = [i for i in xrange(0, len(csum)) if csum[i]!=0]
		#print nonzero_idx

		if len(nonzero_idx) == len(cp.aas02): # all alphabet present in the curren column, drop the gap column (last column)
			nonzero_idx = nonzero_idx[0:len(cp.aas02)-1]
		output_npcmatrix = npcmatrix[:,nonzero_idx]
		aatrack =  ','.join([cp.aas02[j] for j in nonzero_idx])

		#print aatrack
		#print output_npcmatrix
		outfile = '%s_%s.csv' % (outprefix, cols[c])
		np.savetxt(outfile, output_npcmatrix, fmt='%d', delimiter=',')

	cp._info('save %d .csv files' % pfm.msalen)

	'''
	def msacol(self, i):
		return [s[1][i] for s in self.msalist
	'''

if __name__=='__main__':
	cp.dispatch(__name__)

