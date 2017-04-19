import sys
import numpy as np

import common.commp as cp
from protein import protein
import os

# wrapper for multi-processing
def mp_rcg(arglist):
	return rcg(*arglist)

# contact extraction procedure
def rcg(plist, method, cgsize, cutoff, seqdist, title):
	# init stat table
	sm = {}
	AA = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
	for i in xrange(0, len(AA)):
		for j in xrange(i, len(AA)):
			sm['%s%s' % (AA[i], AA[j])] = 0
	if method == 'nearest':
		outfile = '%s.%s.%d.cg' % (title, method, cgsize)
		fout = open(outfile, 'w')
		for p in plist:
			cglist = [] 
			#for ncg in p.contactbynearest(p.atomsbytip('AAtips.py'),cgsize):
			for ncg in p.contactbynearest(p.atoms,cgsize):
				# does not need to sort
				cglist.append('-'.join(['%s%d' % (a.chainID, a.resSeq) for a in ncg]))
				# need sort to match the key in sm
				key = ''.join(sorted([cp.aa2a[a.resName] for a in ncg]))
				sm[key]+=1
			fout.write('%s %s\n' % (p.pdbfile, ' '.join(cglist)))
		fout.close()	

	elif method == 'cutoff': # just for pairwise for now
		outfile = '%s.%s.%d_%d.cg' % (title, method, cutoff, seqdist)
		fout = open(outfile, 'w')
		for p in plist:
			cglist = [] 
			for pcg in p.contactbycutoff(p.atoms,cutoff,seqdist):
				# does not need to sort
				cglist.append('-'.join(['%s%d' % (a.chainID, a.resSeq) for a in pcg]))
				# need sort to match the key in sm
				key = ''.join(sorted([cp.aa2a[a.resName] for a in pcg]))
				sm[key]+=1
			fout.write('%s %s\n' % (p.pdbfile, ' '.join(cglist)))
		fout.close()		

	print os.getpid(), outfile, sum(sm.itervalues())

	return outfile, sm

def main():
	if len(sys.argv) < 3:
		print 'Usage: python proc_rcgstat.py pdblistlist title'
		return

	pdblistlist = sys.argv[1]
	outfile = sys.argv[2]+'.cg.stat'

	#with open(pdblistfile) as fp:
	#	plist = [protein(line.strip()) for line in fp]

	#rcg(plist,'cutoff', 2, 6.5, 6, title)
	arglist=[]
	with open(pdblistlist) as fp:
		for pdblist in fp:
			pdblist = pdblist.strip()

			with open(pdblist) as ffp:
				plist = [protein(line.strip()) for line in ffp]
			arglist.append((plist, 'cutoff', 2, 6.5, 6, pdblist))
			arglist.append((plist, 'nearest', 2, 6.5, 6, pdblist))


	nmp  = len(arglist) if len(arglist) < 20 else 20
	print 'split %d tasks in %d processes\n' % (len(arglist), nmp)

	from multiprocessing import Pool
	pool = Pool(nmp)
	mpret = pool.map(mp_rcg, arglist)

	fout = open(outfile, 'w')
	# each sm is accumulative for each pdblist
	for name, sm in mpret:
		fout.write('%s,%s\n' % (name,' '.join([str(sm[k]) for k in sorted(sm.keys())])))
	fout.write('%s,%s\n' % ('title', repr([k for k in sorted(sm.keys())])))
	fout.close()
	print '\nsave to %s.' % outfile


if __name__ == '__main__':
	main()