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
	if len(sys.argv)<3:
		print 'Usage: python proc_fastaflat.py uniprot.fasta 20905476'
		print 'output: uniprot.fasta.flat'
		exit()

	for p in sys.argv:
		print p,
	print ''

	dbname = sys.argv[1]
	maxnum = int(sys.argv[2])

	fi = fasta_iter(dbname)
	count=0
	fout = open(('%s.flat' % dbname), 'w')
	while(True):
		(head, seq) = fi.next()
		if len(seq)<100:
			continue
		headArray = head.split(' ')
		count+=1
		if count%10000 == 0:
			print '.',
			sys.stdout.flush()
		fout.write('%d %s %s\n' % (len(seq), headArray[0], seq))
		if count == maxnum:
			break
	fout.close()
	print '%d fasta sequences flatened.' % count

if __name__ == '__main__':
	main()