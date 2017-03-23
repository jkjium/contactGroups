import sys
import subprocess as sp 
from collections import defaultdict

# class key put in exclude[] to make class jump
def ck(cathtuple, fold):
	if fold == 1: # jump in C level 
		return '%s' % (cathpuple[1])
	elif fold == 2: # jump in A level
		return '%s %s' % (cathtuple[1], cathtuple[2])
	elif fold == 3: # jump in T level
		return '%s %s %s' % (cathtuple[1], cathtuple[2], cathtuple[3])
	elif fold == 4: # jump in H level
		return '%s %s %s %s' % (cathtuple[1], cathtuple[2], cathtuple[3], cathtuple[4])
	else:
		print 'error: invalid fold-switch-level: %d' % fold
		exit(1)


def main():
	if len(sys.argv)<6:
		print 'Usage: python proc_cathpool.py cath.clf identity% pairs/fold total-num fold-switch-level'
		print 'Usage: python proc_cathpool.py cath.clf 30 20 2000 1'
		return

	dbname = sys.argv[1]
	cutoff = float(sys.argv[2])
	numlimit = int(sys.argv[3])
	total = int(sys.argv[4])
	#the switch fold class level: 1: C, 2:A, 3:T, 4:H'
	fold = int(sys.argv[5])
	outfile = '%s.%dp.%dn.list' % (dbname, cutoff, total)

	for p in sys.argv:
		print p,
	print outfile

	# load cath clf 
	pool = []
	print 'loading seqpool ...'
	sys.stdout.flush()
	with open(dbname) as fp:
		for line in fp: # line format: 1oaiA00 1 10 8 10 1 1 1 1 1 59 1.000
			rline = ' '.join(line.split()) # remove extra space
			strArray = rline.split(' ')
			pool.append((strArray[0], strArray[1], strArray[2], strArray[3], strArray[4]))
	print '%d entries loaded.' % len(pool)

	count=0
	exdict = defaultdict(lambda: 0)
	fout = open(outfile, 'w')

	# main loop
	for i in xrange(0, len(pool)):
		cki = ck(pool[i], 4)
		cke = ck(pool[i], fold)
		for j in xrange(i+1, len(pool)):
			ckj = ck(pool[j], 4)
			if cki != ckj: # make sure structurally similar
				break # forward i

			if exdict[cke]==numlimit:
				print 'swith fold %s' % cke
				break

			# get identity percentile
			ret = sp.check_output(['alnfile.sh',pool[i][0]+'.seq', pool[j][0]+'.seq'])
			pid= float(ret[ret.index('(')+1:ret.index('%')])
			print 'p: %s %s %s %.2f' % (pool[i][0], pool[j][0], cki, pid)
			#print ret

			if pid < cutoff :
				fout.write('%s %s %s %.2f\n' % (pool[i][0], pool[j][0], cki, pid))
				#print 'p: %s %s %s %.2f' % (pool[i][0], pool[j][0], cki, pid)
				count+=1
				exdict[cke]+=1
				if count == total:
					print '%d pairs collected.' % count
					print repr(exdict)
					fout.close()
					exit(0)

	fout.close()

	print '%d pairs collected.' % count
	print repr(exdict)

if __name__ == '__main__':
	main()
