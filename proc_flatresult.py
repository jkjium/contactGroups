import sys
import numpy as np
from alignflat import palign

'''
name = flatArray[0]
program = flatArray[1]
matrix = flatArray[2]
gapopen = float(flatArray[3])
gapextend = float(flatArray[4])
alignlen = int(flatArray[5])
nid = float(flatArray[6])
pid = float(flatArray[7])
nsm = float(flatArray[8])
psm = float(flatArray[9])
ngp = float(flatArray[10])
pgp = float(flatArray[11])
score = float(flatArray[12])
seqAlen = float(flatArray[13])
seqA = flatArray[14]
seqBlen = float(flatArray[15])
seqB = flatArray[16]
'''
# calculate the scores for given flat file
def scoreStr(flatfile, measure):
	nid,pid,nsm,psm,ngp,pgp=[],[],[],[],[],[]	
	with open(flatfile) as fp:
		for line in fp:
			p = palign(line.strip())
			nid.append(p.nid)
			pid.append(p.pid)
			nsm.append(p.nsm)
			psm.append(p.psm)
			ngp.append(p.ngp)
			pgp.append(p.pgp)			

	n = np.array([nid, pid, nsm, psm, ngp, pgp]).T
	print repr(n)
	print repr(np.mean(n, axis=0))
	print repr(np.std(n, axis=0))
	if measure == 'nid':
		ret = np.sum(nid)
	elif measure == 'pid':
		ret = np.sum(pid)
	elif measure == 'nsm':
		ret = np.sum(nsm)
	elif measure == 'psm':
		ret = np.sum(psm)
	elif measure == 'ngp':
		ret = np.sum(ngp)
	elif measure == 'pgp':
		ret = np.sum(pgp)
	elif measure == 'nidm':
		ret = np.mean(nid)
	elif measure == 'pidm':
		ret = np.mean(pid)
	elif measure == 'nsmm':
		ret = np.mean(nsm)
	elif measure == 'psmm':
		ret = np.mean(psm)
	elif measure == 'ngpm':
		ret = np.mean(ngp)
	elif measure == 'pgpm':
		ret = np.mean(pgp)	
	elif measure == 'nids':
		ret = np.std(nid)
	elif measure == 'pids':
		ret = np.std(pid)
	elif measure == 'nsms':
		ret = np.std(nsm)
	elif measure == 'psms':
		ret = np.std(psm)
	elif measure == 'ngps':
		ret = np.std(ngp)
	elif measure == 'pgps':
		ret = np.std(pgp)

	return ret


def main():
	if len(sys.argv) < 5:
		print 'Usage: python proc_flatresult.py prefix smlist flag measure'
		print '[nid,pid,nsm,psm,ngp,pgp] + {m,s}'
		return

	prefix = sys.argv[1]
	smlistfile = sys.argv[2]
	flag = sys.argv[3]
	measure = sys.argv[4]

	# load all the matrix names
	smlist = []
	with open(smlistfile) as fp:
		smlist = fp.readlines()

	# generate suffix (gapopen, gapextent iteration)
	suffix = []
	if flag == '0':
		suffix.append('10_0.5.align.flat')
	else:
		for i in xrange(2,18,2):
			for j in xrange(1,11,2):
				#print '_%d_%.1f.align.flat' % (i, j)
				suffix.append('%d_%.1f.align.flat' % (i, j-0.5))

	# main loop
	# one gap combination with all matrices makes one row
	outfile = '%s.flatscore' % prefix
	fout = open(outfile, 'w')
	count = 0
	for s in suffix: # forward row
		score = []
		for m in smlist: # forward column (one matrix for one column)
			flatfile = '%s_%s_%s' % (prefix,m.strip(),s)
			score.append(scoreStr(flatfile, measure))
		count+=1
		fout.write(','.join(['%.2f' % k for k in score])+'\n')
	print '%d lines saved to: %s' % (count, outfile)
	print repr(score)
	fout.close()

if __name__ == '__main__':
	main()
