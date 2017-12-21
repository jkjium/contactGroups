import sys
import commp as cp

def main():
	if len(sys.argv) < 3:
		cp._info('Usage: python proc_fasep.py fafile out_prefix')
	fafile = sys.argv[1]
	out_prefix = sys.argv[2]

	count=0
	for head, seq in cp.fasta_iter(fafile):
		outfafile = '%s%04d.fa' % (out_prefix, count)
		with open(outfafile, 'w') as fp:
			fp.write('>%s\n%s' % (head, seq))
		cp._info('save to %s' % outfafile)
		count+=1
	cp._info('%d fa files saved.' % count)


if __name__ == '__main__':
	main()