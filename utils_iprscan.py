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


	#Molecular Function: iron ion binding (GO:0005506)
	#Biological Process: oxygen transport (GO:0015671)
	def getgo(self, gostr):
		return set([g.replace(',', '').strip() +')' for g in gostr.split(')') if len(g)!=0])
		

	def dump(self):
		print '------------------------------------'
		print self.name
		for p in self.profile:
			print repr(' '.join(p))
		print '------------------------------------'

	def goset(self):
		goterms = set()
		for p in self.profile:
			if len(p) ==14:
				goterms.update(set([('%s,%s)\n' % (self.name, g.replace(',', '').strip())) for g in p[13].split(')') if len(g)!=0]))
		return list(goterms)


# write a file with name and goterms
def writegoset(arglist):
	if len(arglist) < 2:
		cp._err('Usage: python utils_iprscan.py writegoset iprfile outfile')

	iprfile = arglist[0]
	outfile = arglist[1]

	ipr = iprscan(iprfile)
	gostr = ''.join(ipr.goset())
	if len(gostr)!= 0:
		with open(outfile,'w') as fp:
			fp.write(gostr)
		cp._info('save goterms to %s' % outfile)
	else:
		cp._info('empty goterms in %s' % iprfile)


def foo(arglist):
	ipr = iprscan(arglist[0])
	ipr.dump()
	goterms = ipr.goset()
	for g in goterms:
		print g


if __name__ == '__main__':

	if len(sys.argv)<2:
		cp._err('Usage: python utils_iprscan.py cmd [args ...]')

	dispatch = {
		'writegoset':	writegoset,
		'foo': foo
	}

	if sys.argv[1] not in dispatch:
		cp._err('invalid cmd: %s' % sys.argv[1])
	else:
		dispatch[sys.argv[1]](sys.argv[2:])