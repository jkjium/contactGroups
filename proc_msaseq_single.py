import sys
from itertools import groupby
"""
calculate sequence length of a pfam msa

"""
def fasta_iter(fastafile):
	fh = open(fastafile)
	faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))	
	for header in faiter:
		header = header.next()[1:].strip()
		seq = "".join(s.strip() for s in faiter.next())
		yield header, seq

def main():
	if len(sys.argv) < 2:
		print 'Usage: python proc_msaseq_single.py PF00008_Full.txt'
		return

	msafile = sys.argv[1]
	outfile = msafile[0:7]+'.single.fa'

	fout = open(outfile, 'w')
	for s in fasta_iter(msafile):
		seq = s[1].replace('.','')
		fastr = '>%s\n%s' % (s[0], seq)
		fout.write(fastr)
		break
	fout.close()
	print 'save to %s' % outfile
	

if __name__ == '__main__':
	main()