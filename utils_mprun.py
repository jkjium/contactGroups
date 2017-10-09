import sys
import subprocess
import multiprocessing as mp

import commp as cp

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

	# custom weighting for each task
	# for chunking functionality

	def dump(self):
		print 'task file: %s (%d tasks)' % (self.taskfile, len(self.tasklist))
		print 'use %d processors' % self.ncpu
		#for i in xrange(0, len(self.tasklist)):
		#	print '%d. %s' % (i, self.tasklist[i])


	# mp constant definition
	# mp_print = 0
	# mp_log = 1
	# mp_checkin = 2
	def listener(self):
		print 'listener started ...'
		fout = open(outfile, 'w')
		tcount = 0
		total = len(self.tasklist)
		tstart = time.time()
		while tcount!=total:
			flag, m = q.get()
			print 'worker: %s' % m
			# task finished 
			if flag == cp.mp_checkin:
				tcount+=1
				if tcount%100 == 0:
					timeUsed = int(time.time() - tstart)
					print 'listener: %d/%d in %d seconds, %f hours left.' % (tcount, total, timeUsed, 1.0*total*timeUsed/(tcount*3600))
			# elif flag == cp.mp_info		
			fout.write('%s' % (m))
			fout.flush()				

		fout.close()
		print 'listener: save log file: %s' % outfile

	# does not need watcher
	def exec_pool(self):
		if self.opt not in ['cmd', 'func']:
			print 'invalid opt: %s. process terminated.' % self.opt
			exit()
		print 'mprun(): executing %s with %d processes ...' % (self.opt, self.ncpu)

		# set watcher
		#self.watcher = self.pool.apply_async(self.listener, (self.q, self.logfile))
		# exec func
		async_result = self.pool.map_async(cp.dcall, self.tasklist) if self.opt == 'func' else self.pool.map_async(cp.drun, self.tasklist)
		# waiting all the process
		self.pool.close()
		self.pool.join()		
		for ret in async_result.get():
			if ret!=None:
				print ret
		print 'mprun(): done.' 


	def worker_checkin(self, ret):
		print repr(ret)
		#self.q.put((cp.mp_checkin, ret))



# 
def smp(arglist):
	#
	if len(arglist) < 2:
		print 'Usage: python utils_mprun.py smp {func, cmd}  batch_file.sh ncpu'
		exit()
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
		print 'Usage: python utils_mprun.py cmd pdbfile [args ...]'
		return

	dispatch = {
		'test':test,
		'smp':smp
	}

	if sys.argv[1] in dispatch:
		dispatch[sys.argv[1]](sys.argv[2:])
	else:
		print 'invalid cmd: %s' % sys.argv[1]

if __name__ == '__main__':
	main()