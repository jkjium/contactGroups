import sys
import common.commp as cp
from alignflat import palign
from protein import protein

# calculate RMSD for the aligned positions
def alnRMSD(pa):
	pdbs = pa.pairnames()

	#print repr(pdbs)
	p1 = protein(pdbs[0]+'.aln.pdb')
	rmap1 = cp.posmap(pa.seqA.upper(), p1.seq.upper())

	p2 = protein(pdbs[1]+'.aln.pdb')
	rmap2 = cp.posmap(pa.seqB.upper(), p2.seq.upper())
	'''
	print repr(rmap1)
	print ''
	print repr(rmap2)
	print ''
	'''
	v = []
	w = []
	for i in pa.alnpos():
		p = rmap1[i]
		q = rmap2[i]
		v.append((p1.ca[p].x, p1.ca[p].y, p1.ca[p].z))
		w.append((p2.ca[q].x, p2.ca[q].y, p2.ca[q].z))
		'''
		print 'p1.ca[%d] :' % p,
		p1.ca[p].dump()
		print 'p2.ca[%d] :' % q,
		p2.ca[q].dump()
		print ''
		'''

	return (cp.rmsd(v,w), len(v))

# output the RMSD for all the pair alignments in the flat file
def main():
	if len(sys.argv) < 2:
		print 'Usage: python proc_flatrmsd.py cathpair88.pool.1.align.flat'
		exit(1)

	flatfile = sys.argv[1]
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
		print rmsd,n
		fout.write('%s %d %.4f\n' % (pa.name, n, rmsd))

	fout.close()

if __name__ == '__main__':
	main()