import sys
import operator
import math
from collections import defaultdict
'''
read .pairfreq list and sum up all the items
'''

def pair_freq_sum():
	if len(sys.argv)<2:
		print 'Usage: python proc_pairfreq_sum.py pairfreq.list'
		return

	listfile = sys.argv[1]
	outfile = listfile+'.allpairfreq'

	psm = defaultdict(lambda:0)
	with open(listfile) as fp:
		for pfile in fp:
			pfile = pfile.strip()
			if len(pfile) < 2:
				continue
			with open(pfile) as fp1:
				for line in fp1:
					strarr = line.strip().split(' ')
					if len(strarr) < 2:
						continue
					psm[strarr[0]]+=int(strarr[1])
			print '%s processed.' % pfile

	fo0 = open(outfile+'.0', 'w')
	fo1 = open(outfile+'.1', 'w')
	fo2 = open(outfile+'.2', 'w')

	print 'save to %s.0' % outfile
	print 'save to %s.1' % outfile
	print 'save to %s.2' % outfile

	for k,v in sorted(psm.items(), key=operator.itemgetter(1), reverse=True):
		outstr = '%s %d %.4f\n' % (k, v, math.log10(v))
		if k[0]=='.' and k[1]=='.':
			continue
		if (k[0]!=k[2]) and (k[1]!=k[3]):
			fo2.write(outstr)
		elif (k[0]==k[2] and (k[1]==k[3])):
			fo0.write(outstr)
		else:
			fo1.write(outstr)

	fo0.close()
	fo1.close()
	fo2.close()



def main():
	pair_freq_sum()

if __name__ == '__main__':
	main()