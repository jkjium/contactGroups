from sdii import sdii
from msa import msa
import numpy as np

def main():

	score = np.loadtxt('test_sdii.txt', delimiter=',')
	print score.shape[0]
	sdii_core = sdii(score)
	s = [1,2]
	print sdii_core.w_entropy(sdii_core.data[:,s].T)
	#print sdii_core.w_entropy([0,1])

if __name__ == '__main__':
	main()