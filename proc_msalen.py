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
		print 'Usage: python proc_msalen.py PF00008_Full.txt'
		return

	msafile = sys.argv[1]
	for s in fasta_iter(msafile):
		seq = s[1].replace('.','')
		print '%s %d' % (msafile[0:7], len(seq))
		break

if __name__ == '__main__':
	main()