import sys

"""
continued from proc_p902seedmap.py
convert p90 indices in .allmsacol to seed indices
PF00008_p90.txt.rseq 34 35 100 69 110 111 17 21 22 126 127
output: pfam_seed.allseedcol

"""
def loadmap(mapname):
	posmap={}
	try:
		with open(mapname) as fp:
			for line in fp:
				strarr = line.strip().split(' ')
				posmap[strarr[0]] = strarr[1]
	except IOError:
		print 'error: %s not found.' % mapname

	return posmap

def main():
	if len(sys.argv)!=2:
		print 'Usage: python proc_p902seedmap.py pfam_p90.4.allmsacol'
		print 'output: pfam_p90.4.allseedcol'
		return

	namearr = sys.argv[1].split('.')
	outfile = '%s.allseedcol' % ('_'.join(namearr[:-1]))

	fout = open(outfile, 'w')
	with open(sys.argv[1]) as fp:
		for line in fp:
			strarr = line.strip().split(' ')
			pfam = strarr[0][0:7]
			# PF00008_p902seed.posmap
			posmap = loadmap('%s_p902seed.posmap' % pfam)
			print 'converting %s with %d positions mapping' % (pfam, len(posmap))
			if len(posmap) == 0:
				continue
			seedidx = [posmap[strarr[i]] for i in xrange(1, len(strarr)) if strarr[i] in posmap]
			fout.write('%s %s\n' % (pfam+'.seed', ' '.join(seedidx)))
	fout.close()
	print 'save to: %s' % outfile


if __name__ == '__main__':
	main()