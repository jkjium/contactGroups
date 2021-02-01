import multiprocessing as mp
import numpy as np
import itertools
import math
import os
import sys
import time
import commp as cp 
from utils_pfammsa import pfammsa

def init():
	if len(sys.argv) < 6:
		cp._err('Usage: python mp_worker.py datafile outfile')

	datafile = sys.arg[1]
	outfile = sys.argv[2]

	pfm= pfammsa(datafile)
	cp._info('load %d sequences' % pfm.msanum)

	# calculating total tasks
	totalnum = cp.ncr(len(pfm.msanum), 2)
	tasks = [(i,j) for i in range(0, pfm.msanum) j in range(i+1, pfm.msanum)]
	cp._info('In total %d/%d tasks.' % (len(tasks), totalnum))

	# split tasks into blocks
	tasklist = []
	n = len(tasks)/20 +1
	for i in xrange(0, len(tasks), n):
		tasklist.append(tasks[i:i+n])
	cp._info('spliting tasks into %d blocks' % len(tasklist))
	return (pfm, tasklist, outfile)

# perform needle and get %id %sm %gap
def seq_dist(score, idx_pair):
	name = '%d %d' % (idx_pair[0], idx_pair[1])
	ret = align_exec(name, cmd, matrix, gapopen, gapextend) # needle, B62, 10, 0.5
	'''
	flat = alignparse(name, ret)
	ap = palign(flat)
	return '%s %.2f %.2f %.2f' % (name, ap.pid, ap.psm, ap.pgp)
	'''
	return ret

def worker(score, tasks, q):
	cp._info('worker started.')
	for s in tasks:
		outMessage = seq_dist(score, s)
		q.put(outMessage)
	q.put('done')

def listener(total, outfile, q):
	cp._info('write to file [%s]' % outfile)
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
			if tcount%100 == 0:
				timeUsed = int(time.time() - tstart)
				cp._info('%d/%d in %d seconds, %f hours left.' % (tcount, total, timeUsed, 1.0*total*timeUsed/(tcount*3600)))
			fout.write('%s' % (m))
			fout.flush()
		if tcount == total:
			break
	fout.close()


def main():
	(score, tasklist, outfile) = init()

	manager = mp.Manager()
	q = manager.Queue()
	pool = mp.Pool(22) # cpu_count = 8, 1 for main thread, 1 for listener, 6 for worker
	watcher = pool.apply_async(listener, (sdii_core.totalTask, outfile, q))

	for t in tasklist:
		pool.apply_async(worker, (sdii_core, t, q))

	pool.close()
	pool.join()

if __name__=='__main__':
	main()


