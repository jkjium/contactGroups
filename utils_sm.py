import sys
import commp as cp
import numpy as np

class smatrix(object):
	# read emboss format matrix
	def __init__(self, smfile):
		scorelist = []
		with open(smfile) as fp:
			lines = fp.readlines()

		for line in lines:
			line = line.strip()
			if len(line) < 20: # at least 20 AA
				continue
			if line[0] == '#':
				continue
			if any(c.isdigit() for c in line):
				scorelist.append(line.split())
			else: # alphabet
				aa = line.split()
		#print aa
		# square matrix 'A:E' = 'E:A'
		self.score=dict( ('%s%s' % (aa[i], aa[j-1]), int(scorelist[i][j])) for i in xrange(0,len(scorelist)) for j in xrange(1, len(aa)))
		self.npscore = np.array(scorelist)
		self.aa = aa
		#print repr(self.score)


	def getscore(self, a1, a2):
		if a1 not in self.aa:
			cp._err('Invalid AA: %s' % a1)
		if a2 not in self.aa:
			cp._err('Invalid AA: %s' % a2)

		return self.score['%s:%s' % (a1,a2)]


	def psubdict(self):
		# first get the minimum and translate to remove all the negative values
		minscore = min(self.score.values())
		trans_sm = sorted([(k, self.score[k]-minscore) for k in self.score], key=lambda x: x[1], reverse=True)
		#cp._info(trans_sm)
		for k,v in trans_sm:
			for k1,v1 in trans_sm:
				psubstr = cp.quad_permu([k[0],k1[0],k[1],k1[1]])
				if cp.quadtype(psubstr) == 't2':
					print (psubstr, v, v1, v*v1)


# for Dr. Jernigan 2017 Sep Grant
# mutant entry file as:
#	 Q04771	207	Q	E	
# $ python utils_sm.py batchscore b62 total.del.list mu.del.b62
def proc_batchscore(arglist):
	if len(sys.argv) < 4:
		cp._err('Usage: python utils_sm.py batchscore sm mutant_entyfile')
	sm = smatrix(sys.argv[2])

	with open(sys.argv[3]) as fp:
		mutant = [line.strip().split() for line in fp.readlines() if len(line) > 1]

	scores = [sm.getscore(m[2], m[3]) for m in mutant]
	#print repr(scores)

	'''
	outfile = sys.argv[4]
	with open(outfile, 'w') as fp:
		fp.write(', '.join([str(s) for s in scores]))
	print 'save to %s' % outfile
	'''
	print '= np.array([%s])' % ', '.join([str(s) for s in scores])


	'''
	for i in xrange(0, len(mutant)):
		print '%s %s %s %s : %2i' % (mutant[i][0], mutant[i][1], mutant[i][2], mutant[i][3], scores[i])
	'''


def psub(arglist):
	if len(sys.argv) < 2:
		cp._err('Usage: python utils_sm.py pusb b62 outfile')

	outfile = sys.argv[1]
	sm = smatrix(sys.argv[2])

	sm.psubdict()
	'''
	with open(outfile, 'w') as fp:
		sm.psubdict()
	'''


def test():
	'''
	sm = smatrix('seed7990_full.sm')
	print '%s:%s %d' % (sys.argv[2], sys.argv[3], sm.getscore(sys.argv[2], sys.argv[3]))
	'''
	pass


# main routine
def main():
	if len(sys.argv)<2:
		print 'Usage: python utils_sm.py cmd [args ...]'
		return

	dispatch = {
		'test':test,
		'batchscore': proc_batchscore,
		'psub': psub
	}

	cmd = sys.argv[1]

	if sys.argv[1] not in dispatch:
		print 'invalid cmd: %s' % sys.argv[1]
		return
	else:
		#dispatch[sys.argv[1]]()
		dispatch[sys.argv[1]](sys.argv[2:])

if __name__ == '__main__':
	main()
