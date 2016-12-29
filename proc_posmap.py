import sys
from protein import protein
from atom import atom

def main():
	if len(sys.argv) < 3:
		print 'Usage: python proc_posmap.py pdbfile msaseqfile flatalignfile'
		exit()

	pdbfile = sys.argv[1]
	msaseqfile = sys.argv[2]
	flatalignfile = sys.argv[3]

	print ' '.join(sys.argv)

	gaps = ['.','-']

	# read pdb sequence
	p = protein(pdbfile)
	rpdbseq = p.resArray
	print 'rpdbseq:\n%s' % ','.join(rpdbseq)

	# read msa sequence
	with open(msaseqfile) as fp:
		for line in fp:
			msaseq = line.strip()
			if len(msaseq) < 1:
				print 'error msa seq: %s' % msaseq
	print 'msaseq:\n%s' % msaseq

	'''
	$1:  file name
	$2:  aligned sequence length
	$3:  identity number
	$4:  identity percentile 
	$5:  similarity number
	$6:  similarity percentile
	$7:  gaps number
	$8:  gaps percentile
	$9:  align score
	$10:  seq A pure length
	$11: aligned seq A
	$12: seq B pure length
	$13: aligned seq B
	'''
	with open(flatalignfile) as fp:
		for line in fp:
			alignseq = line.strip()
			alignArray = alignseq.split(' ')
			if len(alignArray) < 13:
				print 'error align seq: %s' % alignseq
	A1 = alignArray[10]
	A2 = alignArray[12]

	print 'align seq a:\n%s' % A1
	print 'align seq b:\n%s' % A2

	#print '\nmatching rpdb with A1'
	p=0
	q=0
	rpdb2A1 = []
	while(True):
		if p==len(rpdbseq) or q==len(A1):
			break
		if rpdbseq[p][1] in gaps:
			p+=1
			continue
		if A1[q] in gaps:
			q+=1
			continue
		#print 'matching rpdb[%d:%s] - A1[%d:%s]' % (p, rpdbseq[p][1].upper(), q, A1[q].upper())
		if rpdbseq[p][1].upper() == A1[q].upper():
			rpdb2A1.append((rpdbseq[p], q))
			p+=1
			q+=1
		else:
			print 'error matching rpdb[%d:%s] - A1[%d:%s]' % (p, rpdbseq[p][1].upper(), q, A1[q].upper())
			exit()

	#print '\nmatching A2 with MSA'
	p=0
	q=0
	A22MSA = {}
	while(True):
		if p==len(A2) or q==len(msaseq):
			break
		if A2[p] in gaps:
			p+=1
			continue
		if msaseq[q] in gaps:
			q+=1
			continue
		#print 'matching A2[%d:%s] - MSA[%d:%s]' % (p, A2[p].upper(), q, msaseq[q].upper())
		if A2[p].upper() == msaseq[q].upper():
			A22MSA[p]=q
			p+=1
			q+=1
		else:
			print 'error matching A2[%d:%s] - MSA[%d:%s]' % (p, A2[p].upper(), q, msaseq[q].upper())
			exit()

	#rpdb2msa = {}
	outfile = pdbfile[0:-5]+'-'+msaseqfile[0:-4]+'.map'
	print 'output to: %s' % outfile
	fout = open(outfile, 'w')
	for (res,pos) in rpdb2A1:
		#rpdb2msa[res] = A22MSA[pos]
		#print '%s-%d -> %d-%s' % (res,pos,A22MSA[pos], msaseq[A22MSA[pos]])
		if pos not in A22MSA:
			# rpdb has more character than uniprot
			fout.write('%s -1 -\n' % res)
		else:
			fout.write('%s %d %s\n' % (res, A22MSA[pos], msaseq[A22MSA[pos]]))
	fout.close()

if __name__ == '__main__':
	main()