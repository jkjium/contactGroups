import sys
import numpy as np
import commp as cp

def mi2dict(sdiifile):
	ret = {}
	with open(sdiifile) as fp:
		for line in fp:
			if len(line) < 2:
				continue
			kv = line.strip().split(' ')
			ret[kv[0]] = float(kv[1])
	return ret

def calc_MIp(arglist):
	"""
	calculate MIp based on .col file and raw sdii(MI)
	for scsc pipeline 2.0
	APC(a,b) = (MIax * MIbx) / avgMI

	"""
	if len(arglist) < 4:
		print 'Usage: python proc_MIp.py colfile sdiifile outfile'
		print 'example: python proc_mip.py PF03176_p90.txt.col PF03176_p90.txt.all_2_sdii PF03176_p90.mip'
		print 'output: p1 p2 mi apc mip'
		return

	colfile = sys.argv[1]
	sdiifile = sys.argv[2]
	outfile = sys.argv[3]
	#print colfile, sdiifile, outfile

	cols = sorted([int(j) for j in np.loadtxt(colfile, delimiter=',')])
	#print repr(cols), len(cols)
	if len(cols) == 0:
		cp._info('err:empty column file: %s' % colfile)
		return

	sdiidict = mi2dict(sdiifile)
	#print sdiidict
	if len(sdiidict) == 0:
		cp._info('err:empty sdii file: %s' % sdiifile)
		return

	avgMI = sum(sdiidict.itervalues()) / len(sdiidict)
	#print 'average MI: %.8f' % avgMI

	MIax = {}
	for i in xrange(0, len(cols)):
		marginMI = 0.0
		for j in xrange(0, len(cols)):
			if i == j:
				continue
			k = ('%s-%s' % (cols[i], cols[j])) if i < j else ('%s-%s' % (cols[j], cols[i]))
			#print k, sdiidict[k]
			marginMI+= sdiidict[k]
		#print marginMI, len(cols)-1
		MIax[cols[i]] = marginMI / (len(cols)-1)
	
	#print repr(MIax)

	MIp = 0.0
	fout = open(outfile, 'w')
	for p in xrange(0, len(cols)):
		for q in xrange(p+1, len(cols)):
			apc = (MIax[cols[p]] * MIax[cols[q]])/avgMI
			k = '%d-%d' % (cols[p], cols[q])
			MIp = sdiidict[k] - apc
			fout.write('%s %.8f %.8f %.8f\n' % (k, sdiidict[k], apc, MIp))
	fout.close()
	print 'save records to %s' % outfile


def main():
	calc_MIp(sys.argv)
if __name__ == '__main__':
		main()	