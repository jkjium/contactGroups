import sys
import commp as cp
import numpy as np

class smatrix(object):
	# read emboss format matrix
	def __init__(self, smfile):
		self.name = smfile
		self.aa = []
		scorelist = []
		with open(smfile) as fp:
			for line in fp:
				line = line.strip()
				if len(line)== 0 or line[0] == '#':
					continue
				if any(c.isdigit() for c in line): # score line
					smline = line.split()
					scorelist.append([int(i) for i in smline[1:]])
				else: # alphabet
					self.aa = line.split()

		npscore = np.array(scorelist)
		if ''.join(self.aa)!= ''.join(cp.smaa1) and ''.join(self.aa)!= ''.join(cp.smaa2):
			cp._err('invalid sm format:\n %s' % self.aa)
		self.core = npscore[:20,:20]
		self.edge = np.copy(npscore)
		self.score =dict(('%s%s' % (self.aa[i],self.aa[j]), self.core[i][j]) for i in xrange(20) for j in xrange(20))

	#
	def dump(self):
		print '%s, max: %d, min: %d' % (self.name, np.max(self.core), np.min(self.core))
		print cp.smstr(self.edge, self.aa)


# combine two matrices with weight
def combinesm(arglist):
	if len(arglist)< 5:
		cp._err('Usage:python utils_sm.py combinesm b62 0.3 scsc 0.7 outfile')

	smfile1 = arglist[0]
	w1 = float(arglist[1])
	smfile2 = arglist[2]
	w2 = float(arglist[3])
	outfile = arglist[4]

	sm1 = smatrix(smfile1)
	sm2 = smatrix(smfile2)
	outcore = sm1.core*w1 + sm2.core*w2
	with open(outfile, 'w') as fp:
		fp.write(outemboss(outcore))
	cp._info('write %s, min: %d, max: %d' % (outfile, np.min(outcore), np.max(outcore)))


# output emboss sm format with b62 edge
def outemboss(core):
	cp.b62edge[:20, :20] = core
	return cp.smstr(cp.b62edge, cp.smaa1)


def test(arglist):
	if len(arglist) < 1:
		cp._err('not enough arglist')
	smfile = arglist[0]
	sm = smatrix(smfile)
	print cp.smstr(sm.core, cp.smaa2)
	print cp.smstr(sm.edge, cp.smaa2)
	print outemboss(sm.core)
	print len(sm.score)
	print sm.score['AC'], sm.score['CA']
	print sm.score['KE'], sm.score['WW']
	sm.dump()


# main routine
def main():
	if len(sys.argv)<2:
		print 'Usage: python utils_sm.py cmd [args ...]'
		return

	dispatch = {
		'combinesm':combinesm,
		'test':test
	}

	if sys.argv[1] not in dispatch:
		cp._err('invalid cmd: %s' % sys.argv[1])
	dispatch[sys.argv[1]](sys.argv[2:])

if __name__ == '__main__':
	main()
