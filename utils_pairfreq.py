import numpy as np
import commp as cp
#from collections import defaultdict
import collections
import itertools
import math

# caluclate weighted background frequencies
# frequencies are normalized into sum of 1
def bgdists(arglist):
	if(len(arglist)<3):
		cp._err('Usage:python utils_pairfreq.py bgeij scorefile weightfile outfile')

	scorefile = arglist[0]
	weightfile = arglist[1]
	outfile = arglist[2]

	score = np.loadtxt(scorefile, delimiter=',')
	nrow = score.shape[0]
	ncol = score.shape[1]
	weight = np.loadtxt(weightfile) # weight is required

	eij = collections.defaultdict(float)
	for i in xrange(0, nrow):
		aafreq = collections.Counter(score[i])
		for k in aafreq:
			eij[k]+=(aafreq[k]*weight[i])
	'''
	print nrow, ncol
	print eij
	'''
	eij.pop(0, None)
	total = sum(eij.values())
	with open(outfile , 'w') as fout:
		#fout.write('%s\n' % ('\n'.join(['bg %s %.8f' % (cp.scoreaa['aa'][k], eij[k]/total) for k in eij])))
		fout.write('%s\n' % ('\n'.join(['bg %s %.8f' % (k, eij[cp.aascore['aa'][k]]/total) for k in cp.aat01])))
	cp._info('save to %s' % outfile)


# calculate weighted frequency
# frequencies each pair of columns are weighted into sum of 1 (no final normalization)
def _wfreq(data, varset, w):
	X = data[:, varset].T
	meff = np.sum(w)
	'''
	print X.T 
	print
	print 'meff: %.f' % meff
	'''
	H = 0
	#print [set(x) for x in X] 
	#print
	prob = {}
	for classes in itertools.product(*[set(x) for x in X]):
		#print '-------------------------'
		#print 'class:' + str(classes)
		v = reduce(np.logical_and, (predictions == c for predictions, c in zip(X, classes)))
		#print "v:", 
		#print v
		if np.sum(v) == 0:
			continue
		'''
		# the original p
		p = np.mean(v) # should divide effective number, which is the sum of all weights
		#print 'p0: ' + str(p)
		prob[classes] = p
		'''
		# apply weights
		p = sum(v*w)/meff
		# print 'pw: ' + str(p)
		prob[classes] = p

	'''
	print
	print '=================='
	print 'probability:'
	print '\n'.join(['%s: %s' % (str(k), str(prob[k])) for k in prob])
	print 'sum: %.2f' % (sum(prob.values()))
	'''
	return prob


# calcuate weighted (correlated) joint frequency in tuple of positions
# output: .tupledist.freq2
def tupledist(arglist):
	if(len(arglist)< 5):
		cp._err('Usage: python utils_pairfreq.py tupledist scorefile columnfile 0,1 weightfile outfile')
	scorefile = arglist[0]
	colfile = arglist[1]
	colidxlist = [int(i) for i in arglist[2].split(',')]
	weightfile = arglist[3]
	outfile = arglist[4]

	score = np.loadtxt(scorefile, delimiter=',')
	weight = np.loadtxt(weightfile) # weight is required
	# scoreaa['aa']
	fout = open(outfile, 'w')
	for line in cp.loadlines(colfile):
		#cols = [int(c) for c in line.split(' ')]
		sarr = line.split(' ')
		cols = [int(sarr[c]) for c in colidxlist]
		wfreq = _wfreq(score, cols, weight)
		for t in wfreq:
			#k = ' '.join(['%d %s' % (i, cp.scoreaa['aa'][i]) for i in t])
			k = ' '.join(['%s' % (cp.scoreaa['aa'][i]) for i in t])
			fout.write('%s %s %.4f\n' % (' '.join([str(i) for i in cols]), k, wfreq[t]))
	fout.close()
	cp._info('save to %s' % outfile)


# calculate qij for current family 
# prepared for the matrix derivation
# input: .tupledist.freq2 from tupledist()
# ['0', '1', 'A', 'C', '0.1667']
# ['0', '1', 'C', 'A', '0.1667']
# ...
# ['2', '3', 'A', 'C', '0.3333']
# ['2', '3', 'C', 'F', '0.5000']
# output: .qij
# aa accumulative_frequency normalized(to 1)_freuency(probability)
# AA 0.13891112, 0.069449
# AC 0.43060278, 0.215280
def tupleqij(arglist):
	if(len(arglist) < 2):
		cp._err('Usage: python tupleqij tuplefile outfile')

	tuplefile = arglist[0]
	outfile = arglist[1]
	plist = [line.split(' ') for line in cp.loadlines(tuplefile)]
	qij = collections.defaultdict(float)
	for c, group in itertools.groupby(plist, lambda x: (x[0], x[1])):
		# for freqs of each different column tuple
		aafreq1 = collections.defaultdict(float)
		aafreq2 = collections.defaultdict(float)
		# get weighed frequency for the two columns
		'''
		print c
		print '================='
		'''
		for e in group:
			aafreq1[e[2]]+=float(e[4])
			aafreq2[e[3]]+=float(e[4])
		'''
		print 'aafreq1:'
		print aafreq1
		print '-----------------------'
		print 'aafreq2:'
		print aafreq2
		print '-----------------------'
		'''
		# the first column
		aafreq1.pop('.', None)
		k = aafreq1.keys()
		k.sort()
		# off-diagonal terms
		for i in xrange(0, len(k)):
			for j in xrange(i+1, len(k)):
				qij['%s%s' % (k[i], k[j])]+=aafreq1[k[i]]*aafreq1[k[j]]
		# diagonal terms
		for i in xrange(0, len(k)):
			qij['%s%s' % (k[i],k[i])]+=(aafreq1[k[i]]*aafreq1[k[i]] / 2.0)
		'''
		print 'qij of aafreq1:'
		print qij
		print '-----------------------'
		'''
		# the second column
		# off-diagonal terms
		aafreq2.pop('.', None)
		k = aafreq2.keys()
		k.sort()
		for i in xrange(0, len(k)):
			for j in xrange(i+1, len(k)):
				qij['%s%s' % (k[i], k[j])]+=aafreq2[k[i]]*aafreq2[k[j]]
		# off-diagonal terms
		for i in xrange(0, len(k)):
			qij['%s%s' % (k[i],k[i])]+=(aafreq2[k[i]]*aafreq2[k[i]] / 2.0)
		'''
		print 'qij of aafreq1 and 2:'
		print qij
		print '-----------------------'
		'''
	aakey = ['%s%s' % (cp.aas01[i],cp.aas01[j]) for i in xrange(0, len(cp.aas01)) for j in xrange(i, len(cp.aas01))]
	total = sum(qij.values())
	if total == 0:
		cp._err('%s zero total' % tuplefile)
	with open(outfile, 'w') as fout:
		fout.write('%s\n' % ('\n'.join([('fg %s %.8f' % (k, qij[k]/total)) for k in aakey])))
		#fout.write('%s\n' % ('\n'.join([('fg %s %.8f' % (k, qij[k]/total)) for k in qij])))
		#fout.write('%s\n' % ('\n'.join([('%s %.8f' % (k, qij[k])) for k in qij])))
	cp._info('save to %s' % outfile)




# combine single frequency and substitution frequency into sm
def wfreq2sm(arglist):
		if len(arglist) < 2:
			cp._err('Usage: python utils_pairfreq.py wfreq2sm combine.wfreq outsmfile')

		wfreqfile = arglist[0] # file contains denominator (background)
		outprefix = arglist[1]

		qij = collections.defaultdict(float)
		eij = collections.defaultdict(float)
		with open(wfreqfile) as fp:
				for line in fp:
						line = line.strip()
						if len(line)==0:
								continue
						sarr = line.split(' ')
						t = sarr[0] # frequency name {bg, fg}
						k = sarr[1] # amino acid name {a, aa}
						f = float(sarr[2]) # weighed frequency value
						# bg M 3511.29 
						# fg AC 0.45
						if t == 'bg':
								eij[k]+=f
						elif t == 'fg':
								qij[k]+=f

		# convert accumulative frequency into probability
		# background probability
		total_e = sum(eij.values())
		for k in eij:
				eij[k]=eij[k]/total_e
		# substitution probability
		total_q = sum(qij.values())
		for k in qij:
				qij[k]=qij[k]/total_q

		# calculate log-odds ratio
		sm = collections.defaultdict(int)
		aakey = ['%s%s' % (cp.aas01[i],cp.aas01[j]) for i in xrange(0, len(cp.aas01)) for j in xrange(i, len(cp.aas01))]
		for k in aakey:
			A = k[0]
			B = k[1]
			if A==B:
				sm[A+B] = int(round(2*math.log(qij[A+B]/(eij[A]*eij[B]),2)))
			else:
				sm[A+B] = int(round(2*math.log(qij[A+B]/(2*eij[A]*eij[B]),2)))
			sm[B+A] = sm[A+B]
		cp._info('sm: %s min: %d, max: %d' % (outprefix, min(sm.values()), max(sm.values())))

		# output emboss format sm
		#embossfile = outprefix + '.emboss.sm'
		embossfile = outprefix
		npemboss = np.array([[sm[A+B] for B in cp.smaa2] for A in cp.smaa2])
		cp.b62edge[:npemboss.shape[0], :npemboss.shape[1]] = npemboss
		with open(embossfile, 'w') as fp:
			fp.write(cp.smstr(cp.b62edge, cp.smaa1))
		cp._info('save sm to %s' % embossfile)


def single_210x400(arglist):
	if len(arglist) < 2:
		cp._err('Usage: python utils_pairfreq.py single_210x400 pairfreq.1 outfile')
	infile = arglist[0]
	outfile = arglist[1]


	dictall = dict( ('%s%s' % (cp.aas01[i],cp.aas01[j]), defaultdict(lambda:0)) for i in xrange(0,len(cp.aas01)) for j in xrange(i+1, len(cp.aas01)) )
	with open(infile) as fp:
		for line in fp:
			line = line.strip()
			if len(line) == 0:
				continue
			sarr = line.split(' ')
			p4 = sarr[0]
			if '.' in p4 or any(c.islower() for c in p4) or len(set(p4).intersection(cp.abaa))!=0:
				continue
			if p4[0] != p4[2]:
				key = [ p4[0], p4[2] ]
				aa = p4[1]
			elif p4[1]!=p4[3]:
				key = [ p4[1], p4[3] ]
				aa = p4[2]
			key.sort()
			#if key=='PP':
			#	print p4
			dictall[''.join(key)][aa]+=float(sarr[2])

	with open(outfile , 'w') as fout:
		for k in dictall:
			nd = cp.rank01(dictall[k])
			H=-sum([(nd[t] * np.log2(nd[t])) for t in nd])
			std = np.std([nd[p] for p in nd])
			#print k, ' '.join(['%s:%.2f' % (a,nd[a]) if a in nd else '%s:0' % (a) for a in cp.aat01])
			fout.write('%.4f %.4f %s %s\n' % (std, H, k, ' '.join(['%.2f' % (nd[a]) if a in nd else '0.00' for a in cp.aat01])))

if __name__ == '__main__':
	cp.dispatch(__name__)
