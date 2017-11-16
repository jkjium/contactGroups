import sys
import os
import time
import subprocess
import multiprocessing as mp

import commp as cp

# mp class for separated jobs
class mprun(object):
	def __init__(self, cmdfile, opt, ncpu=21):
		self.taskfile = cmdfile
		self.logfile = cmdfile + '.mplog'
		self.opt = opt
		self.ncpu = ncpu
		with open(cmdfile) as fp:
			self.tasklist = [line.strip() for line in fp if len(line)>1]
		self.q = mp.Manager().Queue()
		self.pool = mp.Pool(self.ncpu)


	def dump(self):
		print 'task file: %s (%d tasks)' % (self.taskfile, len(self.tasklist))
		print 'use %d processors' % self.ncpu
		#for i in xrange(0, len(self.tasklist)):
		#	print '%d. %s' % (i, self.tasklist[i])

	# does not need watcher
	def exec_pool(self):
		if self.opt not in ['cmd', 'func']:
			cp._err('invalid opt: %s. process terminated.' % self.opt)
		cp._info('mprun(): executing %s with %d processes ...' % (self.opt, self.ncpu))
		# exec func
		async_result = self.pool.map_async(cp.dcall, self.tasklist) if self.opt == 'func' else self.pool.map_async(cp.drun, self.tasklist)
		# waiting all the process
		self.pool.close()
		self.pool.join()		
		for ret in async_result.get():
			if ret!=None:
				cp._info(ret)
		cp._info('mprun(): done.')


# listener process for qrun
def listener(outfile, mpn, taskn, q):
	cp._info('listener started .. ') 
	mpcount = tcount = 0
	tstart = time.time()
	with open(outfile, 'w') as fp:
		#while tcount<taskn:
		while tcount < taskn:
			m = q.get()
			if m == 'done':
				mpcount+=1
			else:
				tcount+=1
				if tcount%100 == 0:
					timeUsed = int(time.time() - tstart)
					cp._info('%d/%d in %d seconds, %f hours left.' % (tcount, taskn, timeUsed, 1.0*taskn*timeUsed/(tcount*3600)))
				fp.write('%s' % (m))
				fp.flush()	
	cp._info('save output to %s' % outfile)


# worker process for qrun
def worker(func, task, q):
	cp._info('worker started.')
	for t in task:
		q.put(func(t))
	q.put('done')

# task split for qrun
def splittask(tasks, n):
		return [tasks[i:i+n] for i in xrange(0, len(tasks), n)]	

# mprun with queue
def qrun(func, taskset, outfile, ncpu):
	q = mp.Manager().Queue()
	taskn = len(taskset)
	if taskn < ncpu:
		ncpu = taskn
	pool = mp.Pool(ncpu+1)
	n = taskn/ncpu + 1
	tasks = splittask(taskset, n)

	pool.apply_async(listener, (outfile, ncpu, taskn, q))
	for tchunk in tasks:
		pool.apply_async(worker, (func, tchunk, q))
	pool.close() #called when the parallelizable part of your main program is finished.
	pool.join()


###################################################################################### 

def smp(arglist):
	if len(arglist) < 2:
		cp._err('Usage: python utils_mprun.py smp {func, cmd}  batch_file.sh ncpu')

	if len(arglist) > 2:
		ncpu = int(arglist[2])
		# {cmdfile, opt, ncpu=13}
		mpr = mprun(arglist[1], arglist[0], ncpu)
	else:
		mpr = mprun(arglist[1], arglist[0])

	mpr.dump()
	mpr.exec_pool()


def test(arglist):
	#mpr = mprun(arglist[0])
	#mpr.dump()
	pass


def main():
	if len(sys.argv)<2:
		cp._err('Usage: python utils_mprun.py cmd pdbfile [args ...]')

	dispatch = {
		'test':test,
		'smp':smp
	}

	if sys.argv[1] in dispatch:
		dispatch[sys.argv[1]](sys.argv[2:])
	else:
		cp._err('invalid cmd: %s' % sys.argv[1])

if __name__ == '__main__':
	main()