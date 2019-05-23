import sys
import commp as cp

def main():
	if len(sys.argv)!=3:
		cp._err('Usage: proc_split.py cath-S20.fa outprefix')
	fafile = sys.argv[1]
	outprefix = sys.argv[2]

	count = 0
	for header, seq in cp.fasta_iter(fafile):
		sarr = header.split(' ')
		with open('%s.%05d.%s.fa' % (outprefix, count, sarr[1]), 'w') as fout:
			fout.write('>%s\n%s\n' % ( header, seq))
		count+=1
	print 'save %d .fa sequences' % count
if __name__=='__main__':
	main()
