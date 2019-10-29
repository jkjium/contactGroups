import sys
import commp as cp
import numpy as np
import copy

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

	def highest(self):
		return self.core.max()
	#
	def dump(self):
		print self.name
		self.edge[:20,:20] = self.core
		print '%s, max: %d, min: %d' % (self.name, np.max(self.core), np.min(self.core))
		print cp.smstr(self.edge, self.aa)

	# output npcore as a file
	# for dendrogram plotting
	def outnpcore(self, outname):
		np.savetxt(outname, self.core, fmt='%.8f', delimiter=' ')

	# copy from below
	def outemboss(self, outname):
		cp.b62edge[:20, :20] = self.core
		with open(outname, 'w') as fout:
			fout.write(cp.smstr(cp.b62edge, cp.smaa1))

	# return an 1x210 np array of current sm
	def out210vec(self):
		return np.array([self.core[i][j] for i in xrange(20) for j in xrange(i,20)])


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

	def stat(self):
		'''
		return # of positive and negative scores
		'''
		p=0
		n=0
		ps=0
		ns=0
		for i in xrange(0,20):
			for j in xrange(i,20):
				if self.core[i][j] > 0:
					p+=1
					ps+=self.core[i][j]
				elif self.core[i][j] < 0:
					n+=1	
					ns+=self.core[i][j]
		return (p,n,ps,ns)

	def translate(self, t):
		'''
		translate score with amount of t
		'''
		for s in self.score:
			self.score[s]+=t
		self.score2core()

	# update sm core with input np.array 1x210
	def updateby210vec(self, np210vec):
		iforward = 0
		for i in xrange(20):
			for j in xrange(i, 20):
				self.core[i][j] = np210vec[iforward]
				self.core[j][i] = np210vec[iforward]
				iforward+=1


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


# compare two matrices 
# save (S1 - S2) 20x20 as a result
# output np readable txt
def diff(arglist):
	if len(arglist)<3:
		cp._err('Usage:python utils_sm.py smfile1 smfile2 outfile(sm1-sm2)')

	smfile1 = arglist[0]
	smfile2 = arglist[1]
	outfile = arglist[2]

	sm1 = smatrix(smfile1)
	sm2 = smatrix(smfile2)
	outcore = sm1.core - sm2.core

	#np.savetxt(outfile, outcore, fmt='%d')
	with open(outfile, 'w') as fout:
		fout.write(outemboss(outcore))
	cp._info('write diff matrix %s' % (outfile))


def gabreed(arglist):
	if len(arglist) < 5:
		errout = [
		"\n",
		"sm cross breeding function",
		"Usage: python utils_sm.py gabreed smlistfile mutate_distribution_file factor_of_mutate factor_of_cross outsmprefix",
		"smlistfile: file contains all the seed sm. one sm for each line",
		"mutate_distribution_file: file contains the probabilities of 0: unchange, +1: score+1, -1: score-1, +2: score+2, -2: score-2 in correponding order",
		"factor_of_mutate: how many times of seed sm with mutation only",
		"factor_of_cross: how many times of seed sm with crossover only"
		"outsmprefix: the prefix for output sm",
		"example: python utils_sm.py gabreed smlist.txt odds.txt 30 sm20190104"
		]
		cp._err('\n'.join(errout))

	smlistfile = arglist[0]
	distfile = arglist[1]
	nm = int(arglist[2])
	ncm = int(arglist[3])
	outprefix = arglist[4]

	# load sm into sm array
	smlist=[]
	with open(smlistfile) as fp:
		for line in fp:
			line = line.strip()
			if len(line)==0:
				continue
			smlist.append(smatrix(line))

	if len(smlist)!=4:
		cp._err('only 4 candidate sm allowed.')

	# load mutation distribution 0 +1 -1 +2 -2
	with open(distfile) as fp:
		for line in fp:
			line = line.strip()
			if len(line) == 0 or line[0] == '#':
				continue
			odds = [float(i) for i in line.split(' ')]
	print repr(odds)

	# mutate loop
	count = 0
	for sm in smlist:
		for i in xrange(nm):
			# get mutate layer 210vec
			#sm.dump()
			overlay = np.random.choice([0,1,-1,2,-2], 210, p=odds)
			#print repr(overlay)
			outsm = copy.deepcopy(sm)
			outsm.updateby210vec(overlay+outsm.out210vec())
			outsm.name = '%s.mu.%02d.sm' % (outprefix, count)
			outsm.outemboss(outsm.name)
			count+=1
			#outsm.dump()
	cp._info('%d mutation sm generated.' % (count))

	# cross over loop
	# each candidate cross nmc times as primary gene
	count=0
	permulist = [[0,1,2,3], [1,0,2,3], [2,1,0,3], [3,1,2,0]]
	for p in xrange(len(smlist)):
		sm = smlist[p]
		permu = permulist[p]
		outsm = copy.deepcopy(sm)
		for k in xrange(ncm): # generate k cross breeds for one sm
			smindex = np.random.choice(permu, 210, p=[0.85, 0.05, 0.05, 0.05])
			idx = 0
			# 210 loop
			for i in xrange(20):
				for j in xrange(i, 20):
					outsm.core[i][j] = smlist[smindex[idx]].core[i][j]
					outsm.core[j][i] = outsm.core[i][j]
					idx+=1
			outsm.name = '%s.mc.%02d.sm' % (outprefix, count)
			outsm.outemboss(outsm.name)
			count+=1
	cp._info('%d crossbreed sm generated.' % (count))

	cp._info('done.')




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
	print '%s' % ' '.join(['%i' % (sm.core[i][j]) for i in xrange(20) for j in xrange(i,20)])


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

def stat(arglist):
	if len(arglist) < 1:
		cp._err('Usage: python utils_sm.py stat b62\noutput: b62.stat')

	sm = smatrix(arglist[0])
	(p,n,ps,ns) = sm.stat()
	outfile = arglist[0]+'.stat'
	with open(outfile, 'w') as fout:
		fout.write('%d %d %d %d\n' % (p,n,ps,ns))
	cp._info('save to %s' % outfile )

# take the arth mean of a set of sm(s)
# for single/double conditional sm
def smmean(arglist):
	if len(arglist) < 2:
		cp._err('Usage: python utils_sm.py smmean smfiles.list outfile')
	smlistfile = arglist[0]
	outfile = arglist[1]
	smlist = [smatrix(smfile) for smfile in cp.loadlines(smlistfile)]

	total = np.zeros((20,20))
	for s in smlist:
		total+=s.core
	outcore = np.rint(total/len(smlist))
	outembossfromcore(outcore, outfile)
	cp._info('mean of %d sm, save to %s' % (len(smlist), outfile))


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

	(p,n,ps,ns)= sm.stat()
	cp._info('p %d n %d ps %d ns %d min: %d max %d' % (p,n,ps,ns,sm.lowest(), sm.highest()))

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
def outembossfromcore(core, outname):
	cp.b62edge[:20, :20] = core
	with open(outname, 'w') as fout:
		fout.write(cp.smstr(cp.b62edge, cp.smaa1))


def outnpcore(arglist):
	if len(arglist) < 2:
		cp._err('Usage: python utils_sm.py smfile outfile')
	smfile = arglist[0]
	outfile = arglist[1]
	sm = smatrix(smfile)
	sm.outnpcore(outfile)
	cp._info('save core to %s' % outfile)


def savecolvec(arglist):
	if len(arglist) < 2:
		cp._err('Usage:python utils_sm.py savecolvec sm.emboss outfile')

	sm = smatrix(arglist[0])
	outfile = arglist[1]

	outstr = ['%s%s %d' % (cp.aat01[i], cp.aat01[j], sm.score[cp.aat01[i]+cp.aat01[j]]) for i in xrange(0, len(cp.aat01)) for j in xrange(i, len(cp.aat01))]
	with open(outfile, 'w') as fout:
		fout.write('%s\n' % ('\n'.join(outstr)))
	cp._info('save to %s' % outfile)


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
		'diff':			diff,
		'gabreed':		gabreed,
		'interpolate': 	interpolate,
		'outblast': 	outblast,
		'outblast62':	outblast62,
		'outnpcore':	outnpcore,
		'printflatten':	printflatten,
		'scalepn':		scalepn,
		'stat':			stat,
		'savecolvec':	savecolvec,
		'smmean':		smmean,
		'transform':	transform,
		'test':			test
	}

	if sys.argv[1] not in dispatch:
		cp._err('invalid cmd: %s' % sys.argv[1])
	dispatch[sys.argv[1]](sys.argv[2:])

if __name__ == '__main__':
	main()
