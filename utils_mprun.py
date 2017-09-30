import sys
import subprocess
import multiprocessing as mp

import commp as cp

class mprun(object):

	def __init__(self, cmdfile, ncpu=13):
		self.taskfile = cmdfile
		self.logfile = cmdfile + '.mplog'
		with open(cmdfile) as fp:
			self.tasklist = [line.split() for line in fp]
		manager = mp.Manager()
		self.q = manager.Queue()
		self.pool = mp.Pool(ncpu)
		self.watcher = self.pool.apply_async(listener, (self.q, self.logfile))

	# custom weighting for each task
	# for chunking functionality

	def dump(self):
		print 'task file: %s' % self.taskfile
		for i in xrange(0, len(self.tasklist)):
			print '%d. %s' % (i, self.tasklist[i])


	# mp constant definition
	# mp_info = 0
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



def test(arglist):
	mpr = mprun(arglist[0])
	mpr.dump()


def main():
	if len(sys.argv)<2:
		print 'Usage: python utils_pfammsa.py cmd pdbfile [args ...]'
		return

	dispatch = {
		'test':test
	}

	if sys.argv[1] not in dispatch:
		print 'invalid cmd: %s' % sys.argv[1]
		return
	else:
		dispatch[sys.argv[1]](sys.argv[2:])

if __name__ == '__main__':
	main()