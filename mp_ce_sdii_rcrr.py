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
	if len(sys.argv) < 5:
		print 'Usage: python mp_ce_sdii_rcrr.py MSATitle weight.file targetVar order'
		print 'Example 1: python mp_ce_sdii_rcrr.py PF07714_full.fa PF07714_full.fa.weight{na} 3128 3'
		print 'Example 1: python mp_ce_sdii_rcrr.py PF07714_full.fa na all 3'
		return

	scoreFile = sys.argv[1]+'.score'
	#rowIndexFile = sys.argv[1]+'.row'
	#colIndexFile = sys.argv[1]+'.col'
	colIndexFile = sys.argv[1]+'.rcol'
	weightfile = sys.argv[2]
	targetVar = sys.argv[3].lower()
	order = int(sys.argv[4])

	print 'score file: [%s]' % scoreFile
	print 'weight file: [%s]' % weightfile
	#print 'row index file: [%s]' % rowIndexFile
	print 'column index file: [%s]' % colIndexFile
	print 'target var: [%s]' % targetVar
	print 'order: [%d]' % order

	outfile = '%s.%s_%d_sdii' % (sys.argv[1], targetVar, order)
	print 'write to [%s]' % outfile

	# msa init
	score = np.loadtxt(scoreFile, delimiter=',')
	w = np.loadtxt(weightfile)

	#rowIndex = [int(i) for i in np.loadtxt(rowIndexFile, delimiter=',')]
	colIndex = [int(j) for j in np.loadtxt(colIndexFile, delimiter=',')]
	print 'reduced data dimension: %s (, %d)' % (repr(score.shape), len(colIndex))
	if len(colIndex) < order:
		print 'exit for insufficient dataset'
		exit()
	varlist = colIndex

	if (targetVar != 'all') and (int(targetVar) not in varlist):
		print 'The alignment for var %s is not significant. exit.' % targetVar
		return 

	# sdii init
	sdii_core = sdii(score)
	print 'Setting varlist to sdii ...'
	sdii_core.setVarlist(varlist) # set sequence weight
	print 'Setting target variable ...'
	sdii_core.setTarget(targetVar)
	print 'Setting task order ...'
	if weightfile != 'na':
		sdii_core.setWeight(w)
		print 'Setting weight ...'
	else:
		print 'No weight applied.'
	sdii_core.setOrder(order)
	print repr(varlist)

	# tasklist init
	# calculating total tasks
	totalnum = cp.ncr(len(varlist), order)
	tasks = []
	if targetVar == 'all':
		print 'generating tasks for all ...'
		for s in set(itertools.combinations(list(range(len(varlist))), order)):
			tasks.append(list(s))
		print 'In total %d/%d for order %d.' % (len(tasks), totalnum, order)
	else:
		print 'generating tasks for variable %s' % targetVar
		for s in set(itertools.combinations(list(range(len(varlist))), order-1)):
			target_idx = varlist.index(int(targetVar))
			if target_idx not in s:
 				st = list(s)
 				st.append(target_idx)
 				tasks.append(st)
		print 'In total %d/%d for order %d.' % (len(tasks), totalnum, order)

	sdii_core.setTotalTask(len(tasks))
	# split tasks into blocks
	tasklist = []
	n = len(tasks)/20 +1
	for i in xrange(0, len(tasks), n):
		tasklist.append(tasks[i:i+n])
	print 'spliting tasks into %d blocks' % len(tasklist)

	print 'init done.'
	return (sdii_core, tasklist, outfile)


def worker(sdii_core, tasks, q):
	print 'worker : %d started.' % os.getpid()
	alphabet = [str(i) for i in sdii_core.varlist]
	for s in tasks:
		#print 'worker: %d: %s          ' % (os.getpid(), '-'.join([(alphabet[i]) for i in s]))
		ret_sdii = sdii_core.calc_sdii(s)
		#outMessage = '[pid:%d] %s %.15f\n' % (os.getpid(), '-'.join([(alphabet[i]) for i in s]), ret_sdii)
		outMessage = '%s %.15f\n' % (' '.join([(alphabet[i]) for i in s]), ret_sdii)
		q.put(outMessage)
	q.put('done')


def listener(total, outfile, q):
	print 'listener started ...'
	print 'listener: write to file [%s]' % outfile
	fout = open(outfile, 'w')
	count = 0
	tcount = 0
	tstart = time.time()
	while True:
		m = q.get()
		if m == 'done':
			count+=1
			print 'listener: %d processes done.' % count
		else:
			tcount+=1
			if tcount%100 == 0:
				timeUsed = int(time.time() - tstart)
				print 'listener: %d/%d in %d seconds, %f hours left.' % (tcount, total, timeUsed, 1.0*total*timeUsed/(tcount*3600))
			#fout.write('%d %s' % (timeUsed, m))
			fout.write('%s' % (m))
			fout.flush()
		if tcount == total:
			break
	fout.close()


def main():

	(sdii_core, tasklist, outfile) = init()


	manager = mp.Manager()
	q = manager.Queue()
#	pool = mp.Pool(mp.cpu_count()) # cpu_count = 8, 1 for main thread, 1 for listener, 6 for worker
	print 'prepare spawn 20 processes ...' 
	pool = mp.Pool(23) # cpu_count = 8, 1 for main thread, 1 for listener, 6 for worker
	watcher = pool.apply_async(listener, (sdii_core.totalTask, outfile, q))


	for t in tasklist:
		pool.apply_async(worker, (sdii_core, t, q))

	pool.close()
	pool.join()

if __name__=='__main__':
	main()


