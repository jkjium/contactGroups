import multiprocessing as mp
import numpy as np
import itertools
import math
import os
import sys
import time
import commp as cp 
from sdii import sdii
#from msa import msa

# varfile is based on the index of .rcol file (from msareduce)

def init():
	if len(sys.argv) < 7:
		cp._err('Usage: python mp_ce_weight.py scorefile colfile weightfile varfile order outfile')

	scorefile = sys.argv[1]
	colfile = sys.argv[2]
	weightfile = sys.argv[3]
	varfile = sys.argv[4]
	order = int(sys.argv[5])
	outfile = sys.argv[6]

	varlist = [int(j) for j in np.loadtxt(colfile, delimiter=',')]
	if len(varlist) < 2:
		cp._err('err: %s not enough varlist: %d' % len(varlist))	

	# msa init
	score = np.loadtxt(scorefile, delimiter=',')
	# sdii init
	sdii_core = sdii(score)
	sdii_core.setVarlist(varlist) 
	# set sequence weight
	if weightfile != 'na':
		w = np.loadtxt(weightfile)
		sdii_core.setWeight(w)
	cp._info('set weight: %r' % sdii_core.isWeighted)
	sdii_core.setOrder(order)

	# calculating total tasks
	tasks = []
	if varfile != 'na':
		with open(varfile) as fp:
			for line in fp:
				line = line.strip()
				if len(line)==0:
					continue
				tasks.append([int(i) for i in line.split(' ')])
		totalnum = len(tasks)
		if totalnum == 0:
			cp._err('task exit for EMPTY TASK.')
		cp._info('In total %d tasks from varfile for order %d.' % (totalnum, order))
	else:
		# generating total tasks
		totalnum = cp.ncr(len(varlist), order)
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
		outMessage = '%s %.8f\n' % (' '.join([(alphabet[i]) for i in s]), ret_sdii)
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
	(sdii_core, tasklist, outfile) = init()

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


