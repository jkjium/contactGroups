import multiprocessing as mp
import numpy as np
import itertools
import math
import os
import sys
import time
import commp as cp 
from sdii import sdii

# calculate entropies for triplet columns
# $ python mp_entropy_mat3.py t.score.10 na t.mp.out
def init():
	if len(sys.argv) < 4:
		cp._err('Usage: python mp_ce_entropy_mat3.py scorefile weightfile outfile')

	scorefile = sys.argv[1]
	weightfile = sys.argv[2]
	outfile = sys.argv[3]

	score = np.loadtxt(scorefile, delimiter=',')
	ncol = score.shape[1]

	# sdii init
	sdii_core = sdii(score)
	#sdii_core.setVarlist(varlist) 
	# set sequence weight
	if weightfile != 'na':
		w = np.loadtxt(weightfile)
		sdii_core.setWeight(w)
	cp._info('set weight: %r' % sdii_core.isWeighted)
	# fixed for the third order
	order = 3
	sdii_core.setOrder(order)

	# generating tasks of variable triplets
	# ncol = 4, s:
	# (0, 1, 2)
	# (0, 2, 3)
	# (1, 2, 3)
	# (0, 1, 3)
	tasks = []
	totalnum = cp.ncr(ncol, order)
	for s in set(itertools.combinations(list(range(ncol)), order)):
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
	for s in tasks:
		idx = ' '.join(['%d' % i for i in s])
		# s = [1,2,3]
		# spectrum_idx_list = [(1,), (2,), (3,), (1, 2), (1, 3), (2, 3), (1, 2, 3)]
		spectrum_idx_list = list(itertools.chain.from_iterable(itertools.combinations(s, i) for i in range(1,4)))
		entropy_spectrum = map(sdii_core.hashed_entropy, spectrum_idx_list)
		outMessage = '%s %s' % (idx, ' '.join(['%.6f' % e for e in entropy_spectrum]))
		q.put(outMessage)
	q.put('done')


def listener(total, outfile, q):
	cp._info('write to file %s' % outfile)
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
				cp._info('%d/%d in %d seconds.' % (tcount, total, timeUsed ))
			fout.write('%s\n' % (m))
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


