#!/usr/bin/python
import numpy as np
import itertools
import math
import time
import sys

from sdii import sdii
from msa import msa


#alphabet = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
alphabet = []
#alphabet = ['X','Y','Z','U','V','W']
#alphabet = ['X1','X2','X3','X4','X5']

def main():
	global alphabet

	if len(sys.argv) < 3:
		print 'Usage: python proc_coevol_sdii.py msafile msapos order'
		print 'Example: python proc_coevol_sdii.py 1k2p_PF07714_seed.txt 472 2'
		return

	msafile = sys.argv[1]
	target = int(sys.argv[2])
	order = int(sys.argv[3])
	print 'msafile: %s' % msafile
	print 'target var: %d' % target
	print 'order: %d' % order

	outfile = '%s.%d_%d_sdii' % (msafile, target, order)
	print 'write to %s' % outfile

	m = msa(msafile)
	score = np.array(m.msaboard())
	print 'data dimension: %s' % repr(score.shape)

	#m.writeScoreboard('1k2p_PF07714_seed.score')
	#return 

	alphabet = [str(i) for i in xrange(0, score.shape[1])]

	sdii_core = sdii(score)
	fout = open(outfile, 'w')
	t0 = time.time()
	for s in set(itertools.combinations(list(range(len(alphabet))), order)): 
	#	if target in s:
	#		ret_sdii = sdii_core.calc_sdii(list(s))
	#		print '%s: %.15f\n' % ('-'.join([(alphabet[i]) for i in s]), ret_sdii)
	#		fout.write('%s %.15f\n' % ('-'.join([(alphabet[i]) for i in s]), ret_sdii))
		ret_sdii = sdii_core.calc_sdii(list(s))
		print '%s: %.15f\n' % ('-'.join([(alphabet[i]) for i in s]), ret_sdii)
		fout.write('%s %.15f\n' % ('-'.join([(alphabet[i]) for i in s]), ret_sdii))

	fout.close()
	t1 = time.time()
	print 'time used: %d seconds' % (t1-t0)


if __name__=="__main__":
	main()
