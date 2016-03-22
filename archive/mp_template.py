import multiprocessing as mp
import os

def worker(arg, q):
	outstr = 'writing by %d - %d' % (arg, os.getpid())
	q.put(outstr)
	q.put('done')


def listener(q):
	count = 0
	while True:
		m = q.get()
		if m == 'done':
			count+=1
		else:
			print 'listener():get: [%s]' % m
		if count == 20:
			print '%d processes done.' % count
			break

def main():
	manager = mp.Manager()
	q = manager.Queue()
	pool = mp.Pool(mp.cpu_count()-2)
	watcher = pool.apply_async(listener, (q,))

	for i in range(20):
		pool.apply_async(worker, (i,q))

	pool.close()
	pool.join()

if __name__=='__main__':
	main()