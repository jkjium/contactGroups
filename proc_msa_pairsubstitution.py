import sys
import numpy as np
#import collections
from collections import defaultdict
import math
from msa import msa
from protein import protein

'''
calculate distinct pair substitution frequency for each pfam
'''

def mp_count(args):
	return calc_psm(*args)

def calc_psm(p,m):
	#psm = defaultdict(lambda:0)
	pair = [int(c) for c in p.split('-')]
	'''
	print ''
	print pair

	print repr(m[:,pair])
	print m[:,pair].shape
	'''

	nrow = m.shape[0]
	'''
	ps = [list(m[i,pair])+list(m[j,pair]) for i in xrange(0, nrow) for j in xrange(i+1, nrow)]
	ups = [permu(list(m[i,pair])+list(m[j,pair])) for i in xrange(0, nrow) for j in xrange(i+1, nrow)]
	for i in xrange(0, len(ps)):
		print repr(ps[i]) + ' -> '+ repr(ups[i])
	'''	
	psm = defaultdict(lambda:0)
	for i in xrange(0, nrow):
		for j in xrange(i+1, nrow):
			psm[permu(list(m[i,pair])+list(m[j,pair]))]+=1

	return psm

# main routine
def pair_substitution():
	if len(sys.argv) < 2:
		print 'Usage: python proc_msa_pairsubstitution.py PF00008_full.txt.rcolpair'
		return

	pairfile = sys.argv[1]
	outfile = pairfile+'.pairfreq'

	psm = defaultdict(lambda:0)
	with open(pairfile) as pf:
		for line in pf: # for each pfam
			strarr = line.strip().split(' ')
			#print strarr
			if len(strarr) < 2:
				print 'skip %s for no top ncg' % pairfile
				continue
			m = msa(strarr[0])
			mm = np.array([list(s[1]) for s in m.msaArray])
			#print mm.shape

			psmlist = []
			for p in strarr[1:]:
				print '%s %s' % (pairfile, p)
				psmlist.append(calc_psm(p,mm))

			for sm in psmlist:
				for k in sm:
					psm[k]+=sm[k]

	# write frequency
	fout = open(outfile,'w')
	for k in psm:
		fout.write('%s %d\n' % (k, psm[k]))
	fout.close()
	print 'save to %s' % outfile




def permu(pairlist):
	'''
	rank = 0  rank = 1  rank = 2  rank = 3
	A 0  C 1  C 0  A 1  D 0  G 1  G 0  D 1
	D 2  G 3  G 2  D 3  A 2  C 3  C 2  A 3
	'''
	rank = pairlist.index(min(pairlist))

	if rank == 1:
		pairlist[0],pairlist[1]=pairlist[1],pairlist[0]
		pairlist[2],pairlist[3]=pairlist[3],pairlist[2]
	elif rank == 2:
		pairlist[0],pairlist[2] = pairlist[2],pairlist[0]
		pairlist[1],pairlist[3] = pairlist[3],pairlist[1]
	elif rank == 3:
		pairlist[0],pairlist[3] = pairlist[3],pairlist[0]
		pairlist[2],pairlist[1] = pairlist[1],pairlist[2]

	#return ''.join(pairlist), rank
	return ''.join(pairlist)


def test():
	# permutation
	'''
	print unify_permutation(['A','C','D','G'])
	print unify_permutation(['C','A','G','D'])
	print unify_permutation(['D','G','A','C'])
	print unify_permutation(['G','D','C','A'])	
	'''
	pass

def main():
	pair_substitution()

if __name__ == '__main__':
	main()
