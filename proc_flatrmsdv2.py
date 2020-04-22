import sys
import commp as cp
from alignflat import palign
from protein import protein

# calculate RMSD for the aligned positions
def alnRMSD(pa):
	pdbs = pa.pairnames()

	#print repr(pdbs)
	p1 = protein('%s.%s.aln.pdb' % (pa.name, pdbs[0]))
	rmap1 = cp.posmap(pa.seqA.upper(), p1.seq.upper())

	p2 = protein('%s.%s.aln.pdb' % (pa.name, pdbs[1]))
	rmap2 = cp.posmap(pa.seqB.upper(), p2.seq.upper())
	
	'''
	print repr(rmap1)
	print ''
	print repr(rmap2)
	print ''
	'''

	v = []
	w = []
	apos = pa.alnpos()
	#print repr(apos)

	if len(apos) == 0 or len(rmap1)==0 or len(rmap2)==0:
		return (0, 0)

	if len(p1.ca)==len(p1.resDict) and len(p2.ca)==len(p2.resDict):
		r1 = p1.ca
		r2 = p2.ca
	else:
		r1 = p1.atomsbygmcenter()
		r2 = p2.atomsbygmcenter()

	#for i in pa.alnpos():
	for i in apos:
		p = rmap1[i]
		q = rmap2[i]
		'''
		v.append((p1.ca[p].x, p1.ca[p].y, p1.ca[p].z))
		w.append((p2.ca[q].x, p2.ca[q].y, p2.ca[q].z))
		'''
		v.append((r1[p].x, r1[p].y, r1[p].z))
		w.append((r2[q].x, r2[q].y, r2[q].z))

		'''
		print 'p1.ca[%d] :' % p,
		p1.ca[p].dump()
		print 'p2.ca[%d] :' % q,
		p2.ca[q].dump()
		print ''
		'''

	return (cp.rmsd(v,w), len(v))

# output the RMSD for all the pair alignments in the flat file
def flatrmsd(arglist):
	if len(arglist) < 1:
		cp._err('Usage: python proc_flatrmsd.py flatrmsd cathpair88.pool.1.align.flat')

	flatfile = arglist[0]
	outfile = '%s.rmsd' % flatfile

	gap = ['.', '-']
	fout = open(outfile, 'w')

	with open(flatfile) as fp:
		lines = fp.readlines()

	fout = open(outfile, 'w')
	for i in xrange(0, len(lines)):
		pa =  palign(lines[i].strip())
		#pa.dump()
		rmsd,n = alnRMSD(pa) 
		#print rmsd,n
		fout.write('%s %d %.4f\n' % (pa.name, n, rmsd))
	cp._info('save to %s' % outfile)
	fout.close()

if __name__ == '__main__':
	cp.dispatch(__name__)