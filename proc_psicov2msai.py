import sys
import commp as cp
'''
append msai to a psicov output file
psicov index: ungapped sequence from its input MSA .aln file (first sequence)
output: msai, shifted aln.index(start from 0), psicov scores
'''
def main():
	if len(sys.argv)!=5:
		cp._err('Usage: proc_psicov2msai.py psicovfile psicov.aln.first.row.seq PF00870_MSA.seq outfile.msai.psicov')
	psicovfile = sys.argv[1]
	alnseqfile = sys.argv[2]
	msaseqfile = sys.argv[3]
	outfile = sys.argv[4]

	alnseq = cp.loadlines(alnseqfile)[0]
	msaseq = cp.loadlines(msaseqfile)[0]

	#rmap[aln.index] = msai
	rmap = cp.posmap_subseq_d(alnseq, msaseq)
	# out format:
	# msai1 msai2 aln.idx1 aln.idx2 psicov.value
	outlist = [] 
	for line in cp.loadlines(psicovfile):
		# aln.1 aln.2 0 8 psicov.value
		# 104   117   0 8 0.678857
		sarr = line.split(' ')
		alnid1= int(sarr[0])-1
		alnid2= int(sarr[1])-1
		msai1 = rmap[alnid1]
		msai2 = rmap[alnid2]
		if alnseq[alnid1]!=msaseq[msai1] or alnseq[alnid2]!=msaseq[msai2]:
			cp._err('rmap error: alnid1: %d %s; msai1: %d %s; alnid2: %d %s; msai2: %d %s' % (alnid1, alnseq[alnid1], msai1, msaseq[msai1], alnid2, alnseq[alnid2], msai2, msaseq[msai2]))
		value = sarr[4]
		outlist.append('%d %d %s %s %s' % (msai1, msai2, alnid1, alnid2, value))

	cp._info('%d pairs mappped. ' % len(outlist))
	with open(outfile,'w') as fout:
		fout.write('%s\n' % ('\n'.join(outlist)))

if __name__=='__main__':
	main()
