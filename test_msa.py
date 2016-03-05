'''
test msa class 
'''
import numpy as np
from protein import protein
from msa import msa

def main():
	#msafile = 'PF07714_seed.txt'
	# test pos map
	'''
	target = '1k2p'
	msafile = '1k2p_PF07714_seed.txt'
	m = msa(msafile, target)
	m.dump()

	p = protein('1k2p.pdb')
	seqi2msai, msai2seqi = m.getPosMap(p)

	s = m.msaArray[0]
	seq = s[1]

	pdbseq = p.seq
	print ''.join([seq[seqi2msai[k]] for k in seqi2msai])
	print ''.join([pdbseq[msai2seqi[l]]for l in msai2seqi])
	'''

	# test msaboard
	msafile = 'test_msa.txt'
	target = '1k2p'
	m = msa(msafile, target)
	msaboard, varlist = m.msaboard(0.0, 0.5)
	print varlist
	print m.seqlen
	print msaboard
	print msaboard.shape
	print m.weight

if __name__ == '__main__':
	main()
