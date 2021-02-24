import multiprocessing as mp
import numpy as np
import itertools
import math
import os
import sys
import time
import commp as cp 
import subprocess as sp 
import utils_testcase as ut
from utils_pfammsa import pfammsa

def init():
	if len(sys.argv) < 3:
		cp._err('Usage: python mp_worker.py datafile outfile')

	datafile = sys.argv[1]
	outfile = sys.argv[2]

	pfm= pfammsa(datafile)
	cp._info('load %d sequences' % pfm.msanum)

	# calculating total tasks
	totalnum = cp.ncr(pfm.msanum, 2)
	tasks = [(i,j) for i in xrange(0, pfm.msanum) for j in xrange(i+1, pfm.msanum)]
	cp._info('In total %d/%d tasks.' % (len(tasks), totalnum))

	# split tasks into blocks
	tasklist = []
	n = len(tasks)/20 +1
	for i in xrange(0, len(tasks), n):
		tasklist.append(tasks[i:i+n])
	cp._info('spliting tasks into %d blocks' % len(tasklist))
	return (pfm, tasklist, outfile)


def align_exec(seqpair, cmd='needle', matrix='B62', gapopen='10.0', gapextend='0.5'):
	#$ ./align.sh needle "ALIGN" "LINE" B62 8 2
	return sp.Popen(['align.sh', cmd, seqpair[0], seqpair[1], matrix, gapopen, gapextend], stdout=sp.PIPE).communicate()[0]
	#return sp.check_output(['align.sh', cmd, seqArr[0], seqArr[1], matrix, gapopen, gapextend])

# perform needle and get %id %sm %gap
def seq_dist(score, idx_pair):
	i = idx_pair[0]
	j = idx_pair[1]
	name = '%d %d' % (i,j)
	seqpair = (score.msalist[i][1], score.msalist[j][1])
	ret = align_exec(seqpair) # needle, B62, 10, 0.5
	flat = ut.alignparse('%d-%d', ret)
	ap = ut.palign(flat)
	return '%s %.4f %.4f %.4f' % (name, ap.pid/100.0, ap.psm/100.0, ap.pgp/100.0)

def worker(score, tasks, q):
	cp._info('worker started.')
	for s in tasks:
		outMessage = seq_dist(score, s)
		q.put(outMessage)
	q.put('done')

def listener(total, outfile, q):
	cp._info('listener::write to file [%s]' % outfile)
	fout = open(outfile, 'w')
	count = 0
	tcount = 0
	tstart = time.time()
	while True:
		m = q.get()
		if m == 'done':
			count+=1
			cp._info('%d processes done.' % count)
		else:
			tcount+=1
			if tcount%1000 == 0:
				timeUsed = int(time.time() - tstart)
				cp._info('%d/%d in %d seconds' % (tcount, total, timeUsed))
			fout.write('%s\n' % (m))
			fout.flush()
		if tcount == total:
			break
	fout.close()


def main():
	(score, tasklist, outfile) = init()

	manager = mp.Manager()
	q = manager.Queue()
	pool = mp.Pool(22) # cpu_count = 8, 1 for main thread, 1 for listener, 6 for worker
	watcher = pool.apply_async(listener, (cp.ncr(score.msanum, 2), outfile, q))

	for t in tasklist:
		pool.apply_async(worker, (score, t, q))

	pool.close()
	pool.join()

if __name__=='__main__':
	main()


