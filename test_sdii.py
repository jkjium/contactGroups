from sdii import sdii
from msa import msa
import numpy as np

def main():

	# test msa weight	
	'''
	msafile = 'test_msa.txt'
	target = '1k2p'
	m = msa(msafile, target)
	score, varlist = m.msaboard(0.0, 0.5)
	print score
	sdii_core = sdii(score)
	print sdii_core.w_entropy(sdii_core.data[:,[0,1]].T)
	sdii_core.setWeight(m.weight)
	print sdii_core.w_entropy(sdii_core.data[:,[0,1]].T)
	print sdii_core.weight
	print 'sum(weight): %f' % sum(sdii_core.weight)
	'''


if __name__ == '__main__':
	main()