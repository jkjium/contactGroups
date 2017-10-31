import sys
import commp as cp
import collections


# iterate .allaafreq file to 
def main():
	if len(sys.argv) < 3:
		cp._err('Usage: python proc_sumaafreq.py 4-scol.allaafreq outfile')

	allaafreqfile = sys.argv[1]
	outfile = sys.argv[2]
	afdict = collections.defaultdict(int)
	total = 0
	cp._info('loading aafreq ...')
	with open(allaafreqfile) as fp:
		# A 9534 0.35175620,V 2712 0.10005903
		for line in fp:
			line = line.strip()
			if len(line) == 0:
				continue
			for f in line.split(','):
				sarr = f.split(' ')
				afdict[sarr[0]]+=int(sarr[1])
				total+=float(sarr[1])

	cp._info('writing %s' % outfile)
	with open(outfile, 'w') as fp:
		for k in afdict:
			fp.write('%s %d %.8f\n' % (k, afdict[k], afdict[k]/total))

	cp._info('save pair substitution: %s' % outfile)

if __name__ == '__main__':
	main()