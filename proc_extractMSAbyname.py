import sys
from itertools import groupby

def fasta_iter(fastafile):
	fh = open(fastafile)
	faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))	
	for header in faiter:
		header = header.next()[1:].strip()
		seq = "".join(s.strip() for s in faiter.next())
		yield header, seq	

def main():
	if len(sys.argv) < 4:
		print 'Usage: python proc_extractMSAbyname.py PF00008_full.txt p90_PF00008.title target PF00008_p90.txt'
		return

	msafile = sys.argv[1]
	titlefile = sys.argv[2]
	target = sys.argv[3]
	outfile = sys.argv[4]

	titles = set()
	with open(titlefile) as fp:
		for line in fp:
			strarr = line.split(' ')
			#print '[%s]' % strarr[0]
			titles.add(strarr[0][1:])

	count=0
	fout = open(outfile ,'w')
	for s in fasta_iter(msafile):
		#print 's0: [%s]' % s[0]
		if (s[0] in titles) or (target in s[0]):
			outstr = '>%s\n%s\n' % (s[0], s[1])
			fout.write(outstr)
			count+=1
	fout.close()
	print 'save %d / %d seqs in %s.' % (count, len(titles), outfile)

if __name__ == '__main__':
	main()