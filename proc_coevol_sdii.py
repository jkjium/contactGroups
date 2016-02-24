#!/usr/bin/python
import numpy as np
import itertools
import math
import time
import sys

from sdii import sdii
from msa import msa
from scipy.special import binom

#alphabet = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
alphabet = []
#alphabet = ['X','Y','Z','U','V','W']
#alphabet = ['X1','X2','X3','X4','X5']

def main():
	global alphabet

	if len(sys.argv) < 4:
		print 'Usage: python proc_coevol_sdii.py msafile cutoff msapos order'
		print 'Example: python proc_coevol_sdii.py 1k2p_PF07714_seed.txt 0.6 472 2'
		return

	msafile = sys.argv[1]
	cutoff = float(sys.argv[2])
	target = sys.argv[3].lower()
	order = int(sys.argv[4])

	print 'msafile: [%s]' % msafile
	print 'cutoff: [%f]' % cutoff
	print 'target var: [%s]' % target
	print 'order: [%d]' % order

	outfile = '%s.%s_%d_sdii' % (msafile, target, order)
	print 'write to [%s]' % outfile

	m = msa(msafile)
	print 'original data dimension: (%d, %d)' % (m.seqNum, m.seqlen)
	score, varlist = m.msaboard(cutoff)
	print 'reduced data dimension: %s' % repr(score.shape)

	alphabet = [str(i) for i in varlist]
	#print alphabet
	#m.writeScoreboard('1k2p_PF07714_seed.score')
	if (target != 'all') and (int(target) not in varlist):
		print 'The alignment for var %s is not significant. exit.' % target
		return 

	pk = binom(len(varlist), order)
	print 'total calculations: %d' % pk
	return
	sdii_core = sdii(score)
	fout = open(outfile, 'w')
	t0 = time.time()
	count = 0
	for s in set(itertools.combinations(list(range(len(alphabet))), order)): 
		if (target == 'all') or (alphabet.index(target) in s):
			count+=1
			ret_sdii = sdii_core.calc_sdii(list(s))
			print ('%d/%d: %s          \n' % (count, pk, '-'.join([(alphabet[i]) for i in s])), end='\r')
			fout.write('%s %.15f\n' % ('-'.join([(alphabet[i]) for i in s]), ret_sdii))

	fout.close()
	t1 = time.time()
	print 'time used: %d seconds' % (t1-t0)


if __name__=="__main__":
	main()
