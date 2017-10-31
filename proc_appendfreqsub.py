import sys
import commp as cp

def main():
	if len(sys.argv) < 4:
		cp._err('Usage: python proc_appendfreqsub.py subfile freqfile outfile')

	subfile = sys.argv[1]
	freqfile = sys.argv[2]
	outfile = sys.argv[3]

	cp._info('load freq')
	fqdict = {}
	with open(freqfile) as fp:
		# $ head *.sum
		# . 1408051 0.03876507
		# A 3204158 0.08821372
		for line in fp:
			line = line.strip()
			if len(line)==0:
				continue
			sarr = line.split(' ')
			fqdict[sarr[0]] = float(sarr[2])


	cp._info('writing sub')
	fout = open(outfile, 'w')
	with open(subfile) as fp:
		# $ head 4-tip-scol.allpsub.ps
		# t9 .SSL 725780 58.970012 0.007024 9
		for line in fp:
			line = line.strip()
			if len(line) == 0:
				continue
			sarr = line.split(' ')
			k = sarr[1]
			nfreqsub = float(sarr[2])/(fqdict[k[0]]*fqdict[k[1]]*fqdict[k[2]]*fqdict[k[3]])
			outstr = '%s %s %s %s %s %.3f %s\n' % (sarr[0], sarr[1], sarr[2], sarr[3], sarr[4], nfreqsub, sarr[5])
			fout.write(outstr)
	fout.close()
	cp._info('save to %s' % outfile)



if __name__ == '__main__':
	main()