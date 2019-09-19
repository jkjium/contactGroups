import sys
import numpy as np
import commp as cp

def mi2dict(sdiifile):
	ret = {}
	stub=[]
	listset=set()
	with open(sdiifile) as fp:
		for line in fp:
			if len(line) < 2:
				continue
			kv = line.strip().split(' ')
			ret['%s %s' % (kv[0], kv[1])] = float(kv[2])
			listset.add(int(kv[0]))
			listset.add(int(kv[1]))
			stub.append('%s %s' % (kv[0], kv[1]))
	colslist = list(listset)
	colslist.sort()
	return ret,stub,colslist

def calc_MIp(arglist):
	"""
	calculate MIp based on .rcol file and raw sdii(MI)
	raw sdii format:
	masi1 msai2 sdii
	APC(a,b) = (MIax * MIbx) / avgMI
	output: single column MIp corresponding to the original order in SDII file (for paste to append)
	"""
	if len(arglist) < 3:
		print 'Usage: python proc_MIp.py colfile sdiifile outfile'
		print 'example: python proc_mip2.py PF03176_p90.txt.all_2_sdii PF03176_p90.mip'
		return

	sdiifile = sys.argv[1]
	outfile = sys.argv[2]
	#print colfile, sdiifile, outfile

	sdiidict, stublist, cols = mi2dict(sdiifile)
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
			k = ('%s %s' % (cols[i], cols[j])) if i < j else ('%s %s' % (cols[j], cols[i]))
			#print k, sdiidict[k]
			if k in sdiidict:
				marginMI+= sdiidict[k]
		#print marginMI, len(cols)-1
		MIax[cols[i]] = marginMI / (len(cols)-1)
	
	#print repr(MIax)
	MIpdict={}
	MIp = 0.0
	for p in xrange(0, len(cols)):
		for q in xrange(p+1, len(cols)):
			apc = (MIax[cols[p]] * MIax[cols[q]])/avgMI
			k = '%d %d' % (cols[p], cols[q])
			if k in sdiidict:
				MIp = sdiidict[k] - apc
				MIpdict[k] = MIp
				print '%s %.8f %.8f %.8f\n' % (k, sdiidict[k], apc, MIp)
			#fout.write('%s %.8f %.8f %.8f\n' % (k, sdiidict[k], apc, MIp))

	with open(outfile, 'w') as fout:
		for k in stublist:
			fout.write('%.8f\n' % MIpdict[k])
	print 'save records to %s' % outfile

def main():
	calc_MIp(sys.argv)
if __name__ == '__main__':
		main()	