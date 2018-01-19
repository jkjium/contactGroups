import sys
import commp as cp


class iprscan(object):
	"""
	data structure of iprscan

	"""
	def __init__(self, iprfile):
		self.name = iprfile
		self.profile = []
		with open(iprfile) as fp:
			for line in fp:
				line = line.strip()
				if len(line) == 0:
					continue
				self.profile.append(line.split('\t'))

	def dump(self):
		print '------------------------------------'
		print self.name
		for p in self.profile:
			print repr(' '.join(p))
		print '------------------------------------'

	def goset(self):
		goterms = set()
		for p in self.profile:
			print len(p)
			if len(p) ==14:
				goterms.update(set(p[13].split(',')))
		print repr(goterms)




def foo(arglist):
	ipr = iprscan(arglist[0])
	ipr.dump()
	ipr.goset()


if __name__ == '__main__':

	if len(sys.argv)<2:
		cp._err('Usage: python utils_iprscan.py cmd [args ...]')

	dispatch = {
		'foo': foo
	}

	if sys.argv[1] not in dispatch:
		cp._err('invalid cmd: %s' % sys.argv[1])
	else:
		dispatch[sys.argv[1]](sys.argv[2:])