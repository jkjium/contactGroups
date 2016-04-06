import multiprocessing as mp
import numpy as np
import itertools
import math
import os
import sys
import time
from sdii import sdii
from msa import msa
from scipy.special import binom


def init():
	if len(sys.argv) < 6:
		print 'Usage: python mproc_coevol_sdii.py msafile weightfile cutoff target_seq msapos order'
		print 'Example 1: python mproc_coevol_sdii.py PF07714_full.fa.r50 PF07714_full.fa.r50.weight 0.6 BTK_HUMAN 3128 3'
		print 'Example 2: python mproc_coevol_sdii.py PF07714_full.fa.s62 NA 0.6 BTK_HUMAN 3128 3'
		return

	msafile = sys.argv[1]
	weightfile = sys.argv[2]
	drop_cutoff = float(sys.argv[3]) # for reduce columns
	targetHeader = sys.argv[4]
	target = sys.argv[5].lower()
	order = int(sys.argv[6])

	print 'msafile: [%s]' % msafile
	print 'weightfile: [%s]' % weightfile
	print 'drop_cutoff: [%f]' % drop_cutoff
	print 'target msa header: [%s]' % targetHeader
	print 'target var: [%s]' % target
	print 'order: [%d]' % order

	outfile = '%s.%s_%d_sdii' % (msafile, target, order)
	print 'write to [%s]' % outfile

	# msa init
	m = msa(msafile)
	m.setTarget(targetHeader)
	print 'original data dimension: (%d, %d)' % (m.seqNum, m.seqlen)
	#weight_cutoff = 0.3 # for weighting msa sequence # taken care of in matlab

	score, varlist = m.msaboard(drop_cutoff) #, weight_cutoff) # return a compact score
	print 'reduced data dimension: %s' % repr(score.shape)

	if (target != 'all') and (int(target) not in varlist):
		print 'The alignment for var %s is not significant. exit.' % target
		return 

	# sdii init
	sdii_core = sdii(score)

	print 'Loading weight ...'
	if weightfile.upper() != 'NA':
		pfam_weight = np.loadtxt(weightfile, delimiter=',')
		print 'Weight vector: %s' % repr(pfam_weight.shape)
	else:
		print 'setting weight: %r' % sdii_core.isWeighted

	print 'Applying weight to sdii data ...'
	sdii_core.setWeight(pfam_weight) # set sequence weight
	print 'Seting varlist to sdii ...'
	sdii_core.setVarlist(varlist) # set sequence weight
	print 'Seting target variable ...'
	sdii_core.setTarget(target)
	print 'Seting task order ...'
	sdii_core.setOrder(order)


	# tasklist init
	# calculating total tasks
	tasks = []
	if target == 'all':
		print 'generating tasks for all ...'
		for s in set(itertools.combinations(list(range(len(varlist))), order)):
			tasks.append(list(s))
		print 'In total %d/%d for order %d.' % (len(tasks), binom(len(varlist), order), order)
	else:
		print 'generating tasks for variable %s' % target
		for s in set(itertools.combinations(list(range(len(varlist))), order-1)):
			target_idx = varlist.index(int(target))
			if target_idx not in s:
 				st = list(s)
 				st.append(target_idx)
 				tasks.append(st)
		print 'In total %d/%d for order %d.' % (len(tasks), binom(len(varlist), order), order)

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
		outMessage = '[pid:%d] %s %.15f\n' % (os.getpid(), '-'.join([(alphabet[i]) for i in s]), ret_sdii)
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
			print 'listener: get %d/%d [%s]' % (tcount, total, m.strip('\n'))
			timeUsed = int(time.time() - tstart)
			fout.write('%d %s' % (timeUsed, m))
			fout.flush()
		if count == 20:
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

	if len(tasklist)!=20:
		print 'mismatch task blocks %d vs number of processes 20' % len(tasklist)
		return

	for t in tasklist:
		pool.apply_async(worker, (sdii_core, t, q))

	pool.close()
	pool.join()

if __name__=='__main__':
	main()


