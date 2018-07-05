import sys
from collections import defaultdict

import commp as cp

# generate sh script for pairing cath sequences with the same T level but different H level
#def pairgen_t1h0(arglist):
def pairgen_t1h0(arglist):
	if len(arglist) < 2:
		cp._err('Usage: python utils_cathdb.py pairgen_t1h0 cath-S20.homolog.list out.sh')

	homologlistfile = arglist[0]
	outfile = arglist[1]

	cathlist = []
	with open(homologlistfile) as fp:
		for line in fp:
			line = line.strip()
			if len(line) == 0:
				continue
			# line format: 107lA00 1 10 530 40
			strarr = line.split(' ')
			pdbname = strarr[0]
			T = '%s %s %s' % (strarr[1], strarr[2], strarr[3])
			H = strarr[4]
			cathlist.append((pdbname, T, H))

	# compair all pdb entry
	with open(outfile, 'w') as fout:
		for i in xrange(0, len(cathlist)):
			for j in xrange(i+1, len(cathlist)):
				p1 = cathlist[i]
				p2 = cathlist[j]
				# same T different H
				if (p1[1]==p2[1]) and (p1[2]!=p2[2]):
					fout.write('paste -d " " %s.seq %s.seq > p.%s.%s.t1h0\n' % (p1[0], p2[0], p1[0], p2[0]))
		cp._info('done. run %s to generate pairfiles' % outfile)

'''
def foo(arglist):
	print arglist
'''

if __name__ == '__main__':
	cp.dispatch(__name__)
