import sys

def main():
	if len(sys.argv) < 4:
		print 'extract tsvfile lines with key make by columns'
		print 'column indices are separated by -'
		print 'python proc_colfilter.py mapfile tsvfile 1-2'
		print 'mapfile:'
		print 'PF00001-ACM3_RAT'
		print 'PF00002-PTH1R_HUMAN'
		print '...'
		return

	mapfile = sys.argv[1]
	tsvfile = sys.argv[2]
	idx = sys.argv[3]

	idxArray = idx.split('-')
	'''
	print 'tsv file: %s' % tsvfile
	print 'column index: %d' % index
	print 'loading map from %s ...' % mapfile
	'''
	fmap = open(mapfile, 'r')
	lines = fmap.readlines()
	fmap.close()

	inmap = []
	for line in lines:
		line = line.strip()
		if len(line)<1: 
			continue
		inmap.append(line)

	#print '%d records loaded' % len(inmap)

	#print 'loading flat record from %s ...' % tsvfile
	ftsv = open(tsvfile, 'r')
	lines = ftsv.readlines()
	ftsv.close()

	for line in lines:
		line = line.strip()
		if len(line) <1:
			continue
		flatArray = line.split(' ')
		key = '-'.join([flatArray[int(i)] for i in idxArray])
		if key  in inmap:
			print '%s' % line

if __name__ == '__main__':
	main()