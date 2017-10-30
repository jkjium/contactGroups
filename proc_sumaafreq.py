import sys
import commp as cp
import collections


# iterate .allaafreq file to 
def main():
	if len(sys.argv) < 2:
		cp._err('Usage: python proc_sumaafreq.py 4-scol.allaafreq')

	allaafreqfile = sys.argv[1]
	afdict = collections.defaultdict(int)
	total = 0
	cp._info('loading aafreq ...')
	with open(allaafreqbfile) as fp:
		# ERNY 15096 t2
		for line in fp:
			if len(line) < 1:
				continue
			sarr = line.strip().split(' ')
			afdict[sarr[0]]+=int(sarr[1])
			total+=float(sarr[1])

	outfile = allaafreqfile + '.sum'
	cp._info('writing %s' % outfile)
	with open(outfile, 'w') as fp:
		for k in afdict:
			fp.write('%s %d %.8f\n' % (k, afdict[k], afdict[k]/total))

	cp._info('save pair substitution: %s' % outfile)

if __name__ == '__main__':
	main()