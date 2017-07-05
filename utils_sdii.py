'''
post-process for sdii result
'''
import sys
import numpy as np
from protein import protein
from sdii import sdii


# read key-value sdii result file
# load into a dictionary
def sdii2dict(sdiifile):
	sdiilist = []
	sdiivalue = []
	fin = open(sdiifile, 'r')
	lines = fin.readlines()
	fin.close()

	print 'sdii2dict(): loading %s' % sdiifile
	for line in lines:
		line = line.strip() # remove the last \n
		if len(line) < 2: # skip empty line
			continue
		strArr = line.split(' ')

		sdiivalue.append(float(strArr[1]))
		sdiilist.append((strArr[0], float(strArr[1])))

	print 'sdii2dict(): %s, %f' % (sdiilist[0][0], sdiilist[0][1])
	print 'sdii2dict(): %s, %f' % (sdiilist[1][0], sdiilist[1][1])
	return sdiilist, sdiivalue


def ssummary(valuelist):
	nplist = np.array(valuelist)
	m = np.mean(nplist)
	s = np.std(nplist, ddof=1)
	print '\nssummary(): ------------------------'
	print 'count:\t%d' % len(valuelist)
	print 'min:\t%f' % np.min(nplist)
	print 'median:\t%f' % np.median(nplist)
	print 'max:\t%f' % np.max(nplist)
	print 'mean:\t%f' % m
	print 'std:\t%f' % s
	print 'ssummary(): ------------------------\n'
	return m, s


def tripletstat():
	if len(sys.argv) < 4:
		print 'tripletstat(): python utils_sdii.py tripletstat triplet_sdii mi_sdii rank_method'
		print 'tripletstat(): python utils_sdii.py tripletstat PF00497_full.txt.all_3_sdii PF00497_full.txt.all_2_sdii quantile'
		print 'tripletstat(): output: PF00497_full.txt.all_3_sdii.quantile'
		return
	sdii3file = sys.argv[2]
	mifile = sys.argv[3]
	rankMethod = sys.argv[4]

	print 'tripletstat(): sdii3file: %s' % sdii3file
	print 'tripletstat(): sdii2file: %s' % mifile
	print 'tripletstat(): rank method: %s' % rankMethod
	print

	miList, miValue = sdii2dict(mifile)
	# mi summary
	ssummary(miValue)

	tripletList, tripletValue = sdii2dict(sdii3file)
	# sdii summary
	M, S = ssummary(tripletValue)

	O = M + S
	normValue = []
	for i in xrange(0, len(tripletValue)):
		if tripletValue[i] < O:
			normValue.append(tripletValue[i])

	print 'tripletstat(): adjust count: %d' % len(normValue)
	normList = np.array(normValue)
 	nM = np.mean(normList)
 	nS = np.std(normList)
	print 'tripletstat(): adjusted min: %f' % np.min(normList)
	print 'tripletstat(): adjusted median: %f' % np.median(normList)
	print 'tripletstat(): adjusted max: %f' % np.max(normList)
	print 'tripletstat(): adjusted mean: %f' % nM
	print 'tripletstat(): adjusted std: %f' % nS
	print

	print 'tripletstat(): filtering informative triplets ...'
	cutoff = nM + 4*nS
	infoList = []
	for k, v in tripletList:
		if v >= cutoff:
			infoList.append((k,v))
	print 'tripletstat(): informative count: %d' % len(infoList)
	print infoList[0]
	print infoList[1]
	print '...'
	print infoList[len(infoList)-1]


def getresset():
	if len(sys.argv) < 2:
		print 'getresset(): python utils_sdii.py getresset result_sdii'
		return

	pdbfile = sys.argv[2]
	chain = sys.argv[3]
	p = protein(pdbfile, chain)
	print p.seq


def dumpdii():
	if len(sys.argv) < 4:
		print 'dumpdii(): python utils_sdii.py dumpdii scorefile 0,1,2'
		print 'varlist must contain at least 3 variables'
		return

	scorefile = sys.argv[2]
	varliststr=sys.argv[3]

	varlist = [int(i) for i in varliststr.split(',')]
	if len(varlist)< 3:
		print 'varlist must contain at least 3 variables'
		return

	print 'varlist: %s' % repr(varlist)

	score = np.loadtxt(scorefile, delimiter=',')
	print 'score: %s' % repr(score)

	s = sdii(score)
	d = s.deltaN(varlist)
	#print repr(d)
	for k in s.deltaN(varlist):
		print '%s: %f' % (k, d[k])


def entropy():
	if len(sys.argv) < 4:
		print 'entropy(): python utils_sdii.py entropy scorefile 0,1'
		return

	scorefile = sys.argv[2]
	varliststr = sys.argv[3]	

	varlist = [int(i) for i in varliststr.split(',')]
	print 'varlist: %s' % repr(varlist)

	score = np.loadtxt(scorefile, delimiter=',')
	print 'score: %s' % repr(score)

	s = sdii(score)
	#print 'entropy [%s]: %f' % (repr(varlist), s.entropy(score[varlist, :]))
	#print 'score%s: %s' % (repr(varlist), repr(score[:,varlist].T))
	print 'entropy [%s]: %f' % (repr(varlist), s.entropy(score[:,varlist].T))




def main():
	if len(sys.argv)<3:
		print 'Usage: python utils_sdii.py cmd pdbfile [args ...]'
		return

	dispatch = {
		'getresset': getresset, 'tripletstat': tripletstat, 'entropy':entropy, 'dumpdii':dumpdii
	}

	cmd = sys.argv[1]

	if cmd not in dispatch:
		print 'invalid cmd string: %s' % cmd
	else:
		dispatch[key]()


if __name__ == '__main__':
	main()
