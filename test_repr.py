import time
def main():
	s = [1,2,3,4,5]

	t0 = time.time()
	for i in xrange(0, 10000000):
		s_str = repr(s)
	t1 = time.time()
	print 'repr: %d' % (t1-t0)	


	t0 = time.time()
	for i in xrange(0, 10000000):
		s_str = str(s)
	t1 = time.time()
	print 'str: %d' % (t1-t0)	


if __name__=='__main__':
	main()
