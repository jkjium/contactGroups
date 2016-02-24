'''
test msa class 
'''
import numpy as np
from protein import protein
from msa import msa

def main():
	#msafile = 'PF07714_seed.txt'
	msafile = 'test_msa.txt'
	m = msa(msafile)
	m.dump()

	# test pos map
	'''
	msafile = '1k2p_PF07714_seed.txt'
	p = protein('1k2p.pdb')
	posmap = m.getPosMap(p)
	s = m.msaArray[0]
	seq = s[1]
	print ''.join([seq[posmap[k]] for k in posmap])
	'''

	# test msaboard
	msaboard = np.array(m.msaboard())
	print msaboard
	print msaboard.shape

if __name__ == '__main__':
	main()
