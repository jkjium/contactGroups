import sys
import commp as cp
import numpy as np

class smatrix(object):
	# read emboss format matrix or 20x20 core matrix
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
		if ''.join(self.aa)!= ''.join(cp.smaa1) and ''.join(self.aa[:20])!= ''.join(cp.smaa2):
			cp._err('invalid sm format:\n %s' % self.aa)
		self.core = npscore[:20,:20]
		self.edge = np.copy(npscore)
		self.score =dict(('%s%s' % (self.aa[i],self.aa[j]), self.core[i][j]) for i in xrange(20) for j in xrange(20))

	def lowest(self):
		return self.core.min()
	#
	def dump(self):
		self.edge[:20,:20] = self.core
		print '%s, max: %d, min: %d' % (self.name, np.max(self.core), np.min(self.core))
		print cp.smstr(self.edge, self.aa)

	def score2core(self):
		'''
		update self.score according to new score dictionary
		'''
		scoreupdate = []
		for i in xrange(20):
			scoreupdate.append([self.score['%s%s' %(self.aa[i],self.aa[j])] for j in xrange(20)])
		self.core = np.array(scoreupdate)


	def scale_pn(self, p, n):
		'''
		scaling score by positive and negative scores
		p: scaling for positive score
		n, scaling for negative score
		'''
		for s in self.score:
			i = self.score[s]
			self.score[s] = i*p if i>0 else i*n
		self.score2core()

	def translate(self, t):
		'''
		translate score with amount of t
		'''
		for s in self.score:
			self.score[s]+=t
		self.score2core()




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





def interpolate(arglist):
	if len(arglist) < 3:
		cp._err('Usage: python utils_sm.py interpolate SCSC b62')
	smfile1 = arglist[0]
	smfile2 = arglist[1]

	sm1 = smatrix(smfile1)
	sm2 = smatrix(smfile2)

	neq=[]
	nne=[]
	for i in xrange(20):
		for j in xrange(i, 20):
			if sm1.core[i][j] == sm2.core[i][j]:
				neq.append((i,j))
			else:
				nne.append((i,j))

	print '%d identical, %d different' % (len(neq), len(nne))
	#print [(sm1.aa[i], sm1.aa[j]) for (i,j) in neq]


def printflatten(arglist):
	'''
	print 210 elements into a list
	'''
	if len(arglist) < 1:
		cp._err('Usage: python utils_sm.py printflatten b62')
	smfile = arglist[0]
	sm = smatrix(smfile)
	print '%s\n' % ' '.join(['%i' % (sm.core[i][j]) for i in xrange(20) for j in xrange(i,20)])


def scalepn(arglist):
	if len(arglist) < 4:
		# scale B62 positive score by 2 and leave negative unchanged
		# result save to outsm
		cp._err('Usage: python utils_sm.py scalepn b62 2 1 outsm')

	smfile = arglist[0]
	p = float(arglist[1])
	n = float(arglist[2])
	outsm = arglist[3]

	sm = smatrix(smfile)
	sm.scale_pn(p,n)

	with open(outsm, 'w') as fout:
		fout.write(outemboss(sm.core))
	cp._info('save sm to %s, min: %d, max: %d' % (outsm, np.min(sm.core), np.max(sm.core)))


def transform(arglist):
	'''
	transform sm by 1. translate; 2. scale p,n
	'''
	if len(arglist) < 6:
		cp._err('Usage: python utils_sm.py transform b62 10 2 1 0 outsm')

	smfile = arglist[0]
	t = float(arglist[1])
	p = float(arglist[2])
	n = float(arglist[3])
	order = int(arglist[4])
	outsm = arglist[5]


	cp._info('\ntranslate: %.2f\nscale: + %.2f, - %.2f\norder: %d (0: translate then scale; 1: scale then translate)' % (t,p,n,order))
	sm = smatrix(smfile)
	if order == 0:
		sm.translate(t)
		sm.scale_pn(p,n)
	else:
		sm.scale_pn(p,n)
		sm.translate(t)

	sm.dump()

	with open(outsm, 'w') as fout:
		fout.write(outemboss(sm.core))
	cp._info('save sm to %s, min: %d, max: %d' % (outsm, np.min(sm.core), np.max(sm.core)))



def outblast(arglist):
	if len(arglist) < 1:
		cp._err('Usage: utils_sm.py outblast smemboss')

	scsc = smatrix(arglist[0])
	prefix = '#include <util/tables/raw_scoremat.h>\n/* %s */\nstatic const TNCBIScore s_Blosum80PSM[25 * 25] = {' % arglist[0]
	suffix = '};\nconst SNCBIPackedScoreMatrix NCBISM_Blosum80 = {\n    "ARNDCQEGHILKMFPSTWYVBJZX*",\n    s_Blosum80PSM,\n    -6\n};\n'
	cp.b80blast[:20, :20] = scsc.core
	blastsm = ',\n'.join(['\t%s' % (','.join(['%3i' % n for n in cp.b80blast[i,:]])) for i in xrange(len(cp.aablast))])
	with open('sm_blosum80.c.sub','w') as fp:
		fp.write('\n'.join([prefix, blastsm, suffix]))
	cp._info('save %s to sm_blosum80.c.sub' % (arglist[0]))


def outblast62(arglist):
	if len(arglist) < 1:
		cp._err('Usage: utils_sm.py outblast62 smemboss')

	scsc = smatrix(arglist[0])
	prefix = '#include <util/tables/raw_scoremat.h>\n/* %s */\nstatic const TNCBIScore s_Blosum62PSM[25 * 25] = {' % arglist[0]
	#suffix = '};\nconst SNCBIPackedScoreMatrix NCBISM_Blosum62 = {\n    "ARNDCQEGHILKMFPSTWYVBJZX*",\n    s_Blosum62PSM,\n    -4\n};\n'
	suffix = '};\nconst SNCBIPackedScoreMatrix NCBISM_Blosum62 = {\n    "ARNDCQEGHILKMFPSTWYVBJZX*",\n    s_Blosum62PSM,\n    %d\n};\n' % (scsc.lowest())
	cp.b62blast[:20, :20] = scsc.core
	blastsm = ',\n'.join(['\t%s' % (','.join(['%3i' % n for n in cp.b62blast[i,:]])) for i in xrange(len(cp.aablast))])
	with open('sm_blosum62.c.sub','w') as fp:
		fp.write('\n'.join([prefix, blastsm, suffix]))
	cp._info('save %s to sm_blosum62.c.sub' % (arglist[0]))	


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
		'combinesm':	combinesm,
		'interpolate': 	interpolate,
		'outblast': 	outblast,
		'outblast62':	outblast62,
		'printflatten':	printflatten,
		'scalepn':		scalepn,
		'transform':	transform,
		'test':			test
	}

	if sys.argv[1] not in dispatch:
		cp._err('invalid cmd: %s' % sys.argv[1])
	dispatch[sys.argv[1]](sys.argv[2:])

if __name__ == '__main__':
	main()
