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


# calculate RMSD based on alnpos of pA (protein 1)
def crossRMSD(pA, pB):
	pdbs = pA.pairnames()

	p1 = protein(pdbs[0]+'.aln.pdb')
	m2sA1 = cp.posmap(pA.seqA.upper(), p1.seq.upper()) # input alignment position gives residue index

	pdbA1ResIndex = [m2sA1[i] for i in pA.alnpos()] # get residue index of protein 1 from alignment A

	s2mB1 = cp.posmap(p1.seq.upper(), pB.seqA.upper()) #input residue index gives alignment position

	alnposB = [s2mB1[i] for i in pdbA1ResIndex] # aligned positions in alignment B based on alignment A

	p2 = protein(pdbs[1]+'.aln.pdb')
	m2sB2 = cp.posmap(pB.seqB.upper(), p2.seq.upper()) # given position gives resi index in align B

	for i in alnposB:
		p = pdbA1ResIndex[i]
		q = m2sB2[i]
		v.append((p1.ca[p].x, p1.ca[p].y, p1.ca[p].z))
		w.append((p2.ca[q].x, p2.ca[q].y, p2.ca[q].z))

	return cp.rmsd(v,w)


# output the cross RMSD for all the pair alignments in the flat file
def main():
	if len(sys.argv) < 3:
		print 'Usage: python proc_flatrmsd.py b62.flat cb2.flat output.crmsd'
		exit(1)

	flatA = sys.argv[1]
	flatB = sys.argv[2]
	outfile = sys.argv[3]

	gap = ['.', '-']
	fout = open(outfile, 'w')

	with open(flatA) as fp:
		linesA = fp.readlines()

	with open(flatB) as fp:
		linesB = fp.readlines()

	if len(linesA)!=len(linesB):
		print 'error: flat files line number does not match\n %s: %d - %s: %d' % (flatA, flatB, len(liensA), len(linesB))
		return
	
	fout = open(outfile, 'w')
	for i in xrange(0, len(linesA)):
		pA =  palign(linesA[i].strip())
		pB =  palign(linesB[i].strip())
		#pa.dump()
		rmsdA,nA = alnRMSD(pA) 
		rmsdB,nB = alnRMSD(pB) 

		if pA.name != pB.name:
			print 'error: unmatched name %s - %s' % (pA.name, pB.name)
			return

		rmsdbaseA = crossRMSD(pA, pB)
		rmsdbaseB = crossRMSD(pB, pA)

		fout.write('%s %d %.4f %d %.4f, %.4f %.4f\n' % (pA.name, nA, rmsdA, nB, rmsdB, rmsdbaseA, rmsdbaseB))

	fout.close()

if __name__ == '__main__':
	main()