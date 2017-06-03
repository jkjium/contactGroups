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
	# U!!!
	AA = [
		'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 
		'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'U','V', 'W', 'Y',
		]
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
		outfile = '%s.%s.%s_%d.cg' % (title, method, str(cutoff), seqdist)
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
	if len(sys.argv) < 2:
		print 'Usage: python proc_aastat.py pdblist'
		return

	pdblist = sys.argv[1]
	outfile = sys.argv[1]+'.aa.stat'

	AA = ['G','P','C','U','A','I','L','V','M','F','W','N','Q','S','T','Y','D','E','R','H','K']
	# init scoreboard
	score = {}
	for a in AA:
		score[a]=0

	with open(pdblist) as fp:
		#plist = [protein(line.strip()) for line in fp]
		for line in fp:
			p = protein(line.strip())
			print p.pdbfile
			for i in xrange(0, len(p.seq)):
				score[p.seq[i]]+=1

	fout = open(outfile, 'w')
	fout.write(' '.join([str(score[aa]) for aa in AA]))
	fout.close()
	print 'save to %s.' % outfile

if __name__ == '__main__':
	main()