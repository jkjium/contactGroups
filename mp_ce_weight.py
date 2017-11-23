import multiprocessing as mp
import numpy as np
import itertools
import math
import os
import sys
import time
import commp as cp 
from sdii import sdii
from msa import msa


def init():
	if len(sys.argv) < 6:
		cp._err('Usage: python mp_ce_sdii_rcrr.py scorefile colfile weightfile order outfile')

	scorefile = sys.argv[1]
	colfile = sys.argv[2]
	weightfile = sys.argv[3]
	order = int(sys.argv[4])
	outfile = sys.argv[5]

	# msa init
	score = np.loadtxt(scorefile, delimiter=',')
	varlist = [int(j) for j in np.loadtxt(colfile, delimiter=',')]
	w = np.loadtxt(weightfile)

	if len(varlist) < 2:
		cp._err('err:not enough varlist: %d' % len(varlist))

	# sdii init
	sdii_core = sdii(score)
	sdii_core.setVarlist(varlist) # set sequence weight
	sdii_core.setWeight(w)
	sdii_core.setOrder(order)

	# calculating total tasks
	totalnum = cp.ncr(len(varlist), order)
	tasks = []
	for s in set(itertools.combinations(list(range(len(varlist))), order)):
		tasks.append(list(s))
	cp._info('In total %d/%d for order %d.' % (len(tasks), totalnum, order))

	sdii_core.setTotalTask(len(tasks))
	# split tasks into blocks
	tasklist = []
	n = len(tasks)/20 +1
	for i in xrange(0, len(tasks), n):
		tasklist.append(tasks[i:i+n])
	cp._info('spliting tasks into %d blocks' % len(tasklist))

	return (sdii_core, tasklist, outfile)


def worker(sdii_core, tasks, q):
	cp._info('worker started.')
	alphabet = [str(i) for i in sdii_core.varlist]
	for s in tasks:
		ret_sdii = sdii_core.calc_sdii(s)
		outMessage = '%s,%.8f\n' % (','.join([(alphabet[i]) for i in s]), ret_sdii)
		q.put(outMessage)
	q.put('done')


def listener(total, outfile, q):
	cp._info('listener: write to file [%s]' % outfile)
	fout = open(outfile, 'w')
	count = 0
	tcount = 0
	tstart = time.time()
	while True:
		m = q.get()
		if m == 'done':
			count+=1
			cp._info('listener: %d processes done.' % count)
		else:
			tcount+=1
			if tcount%100 == 0:
				timeUsed = int(time.time() - tstart)
				cp._info('listener: %d/%d in %d seconds, %f hours left.' % (tcount, total, timeUsed, 1.0*total*timeUsed/(tcount*3600)))
			#fout.write('%d %s' % (timeUsed, m))
			fout.write('%s' % (m))
			fout.flush()
		#if count == 20:
		if tcount == total:
			break
	fout.close()


def main():

	(sdii_core, tasklist, outfile) = init()

	manager = mp.Manager()
	q = manager.Queue()
	cp._info('prepare spawn 20 processes ...')
	pool = mp.Pool(22) # cpu_count = 8, 1 for main thread, 1 for listener, 6 for worker
	watcher = pool.apply_async(listener, (sdii_core.totalTask, outfile, q))

	#if len(tasklist)!=20:
	#	cp._err('mismatch task blocks %d vs number of processes 20' % len(tasklist))

	for t in tasklist:
		pool.apply_async(worker, (sdii_core, t, q))

	pool.close()
	pool.join()

if __name__=='__main__':
	main()


