import sys

def main():
	if len(sys.argv) < 4:
		print 'check tsvfile[index] with key in mapfile and append value in tsvfile'
		print 'python proc_addcolumn.py mapfile tsvfile 0_start_column_index'
		return

	mapfile = sys.argv[1]
	tsvfile = sys.argv[2]
	index = int(sys.argv[3])

	'''
	print 'tsv file: %s' % tsvfile
	print 'column index: %d' % index
	print 'loading map from %s ...' % mapfile
	'''
	fmap = open(mapfile, 'r')
	lines = fmap.readlines()
	fmap.close()

	inmap = {}
	for line in lines:
		line = line.strip()
		if len(line)<1: 
			continue
		mapArray = line.split(' ')
		inmap[mapArray[0]] = mapArray[1]

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
		if flatArray[index] not in inmap:
			print 'error: no map for %s' % flatArray[index]
		else:
			print '%s %s' % (line, inmap[flatArray[index]])

if __name__ == '__main__':
	main()