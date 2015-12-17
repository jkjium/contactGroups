#!/usr/bin/python
import numpy as np
import itertools
import math
import csv
import time
import sys

def entropyS(X):
	n_instances = len(X[0])
	#print n_instances
	H = 0
	for classes in itertools.product(*[set(x) for x in X]):
		#print classes
		v = np.array([True] * n_instances)
		for predictions, c in zip(X, classes):
			v = np.logical_and(v, predictions == c)
		p = np.mean(v)
		H += -p * np.log2(p) if p > 0 else 0
	return H

def entropy2(X, Y):
	probs = []
	for c1 in set(X):
		for c2 in set(Y):
			probs.append(np.mean(np.logical_and(X == c1, Y == c2)))
	#print probs
	return np.sum((-p * np.log2(p) if p > 0 else 0) for p in probs)

# tuple of column-wise observations for all the variables
# x1 x2
# 0	 0
# 0  1
# 1  1
# list of all combinations:
# Prob(0, 0) = 1/3
# Prob(0, 1) = 1/3
# Prob(1, 0) = 0
# Prob(1, 1) = 1/3 

# entropy(x1, x2) = 1.585
def entropy(X):
	return np.sum(-p * np.log2(p) if p > 0 else 0 for p in (np.mean(reduce(np.logical_and, (predictions == c for predictions, c in zip(X, classes)))) for classes in itertools.product(*[set(x) for x in X])))

def entropy_single(X):
	probs = [np.mean(X == c) for c in set(X)]
	#print probs
	return np.sum(-p * np.log2(p) for p in probs)

def II(varset, data):
	iiv=0.0
	n = len(varset)
	subsets = list(itertools.chain(*[itertools.combinations(range(n), i) for i in range(n+1)]))
	for s in subsets:
		if len(s) == 0:
			continue
		else:
			iiv+=math.pow(-1, len(s))*entropy(data[:,s].T)
	return -iiv

# dependence calculation
def deltaN_bar(varset, data):
	deltaX_bar = 1.0
	n = len(varset)
	#print 'deltaN_bar(): varset: \n', varset 
	#print 'deltaN_bar(): data: \n' , data

	for index in varset:
		#print '\n\ndeltaN_bar(): index: %d\n' % index
		deltaX = 0.0
		#subsets = list(itertools.chain(*[itertools.combinations(range(n), i) for i in range(n+1)]))
		subsets = list(itertools.chain.from_iterable(itertools.combinations(varset, i) for i in range(len(varset)+1)))
		for tau in subsets:
			if index in tau:
				#print '\ndeltaN_bar(): tau in subsets: ' , tau
				#print 'deltaN_bar(): data[:,tau].T: ' , data[:,tau].T
				#print 'deltaN_bar(): H(tau): %f' % entropy(data[:,tau].T)
				# differential interaction information
				#print 'deltaN_bar(): deltaX: %f' % (math.pow(-1, len(tau)+1)*entropy(data[:,tau].T))
				deltaX+=math.pow(-1, len(tau)+1)*entropy(data[:,tau].T)
				#print 'deltaN_bar(): cumulative deltaX: %f' % deltaX
		deltaX_bar*=deltaX
	return math.pow(-1, n)*deltaX_bar




def main():
	'''
	x0 = np.array([1,1,2,2,3])
	x1 = np.array([0,0,1])
	x2 = np.array([0,1,1])
	x3 = np.array([x1,x2])

	print 'entropy_signle(x1):\t%f' % entropy_single(x1)
	print 'entropy_single(x2):\t%f' % entropy_single(x2)
	print 'entropy([x1]):   \t%f' % entropy([x1])
	print 
	print 'entropy2(x1,x2): \t%f' % entropy2(x1,x2)
	print 'entropyS([x1, x2]): \t%f' % entropyS([x1,x2])
	print 'entropy([x1,x2]): \t%f' % entropy([x1,x2])
	print
	print x3
	print 'entropyS(x3): \t%f' % entropyS(x3)
	print 'entropy(x3): \t%f' % entropy(x3)

	n = 3
	k = 1
	#subsets = list(itertools.chain(*[itertools.combinations(range(n), i) for i range(n+1)]))
	subsets = list(itertools.chain(*[itertools.combinations(range(n), i) for i in range(n+1)]))
	print subsets
	for s in subsets:
		if k in s:
			 print len(s), s

	filename = 'test.txt'
	data = np.loadtxt(filename, delimiter=',')
	print data
	print data[:,0] # get the first column
	tau = [0,1]
	print data[:,tau].T
	print
	print 'start deltaN_bar'

	varset=[0,1]
	print 'deltaN_bar: %f' % deltaN_bar(varset, data)
	'''

	'''
	print 'new subset:'
	subsets = list(itertools.chain.from_iterable(itertools.combinations(varset, i) for i in range(len(varset)+1)))
	for s in subsets:
		print s
	'''
'''
	varset=[0,1]
	filename = 'test.txt'
	data = np.loadtxt(filename, delimiter=',')
	print '\n\ninteraction information: %f' % II(varset, data)
	#print 
	#print
	#print 'Delta K\n'
	# load score data
'''
	scorefile = '5A.hcg.score'
	score = np.loadtxt(scorefile, delimiter=',')
	#print 'score:\n', score
	#print 'AC: %f' % II([0,1], score)

	alphabet = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
	fout = open(scorefile+'.result', 'w')
	t1 = time.time()
	'''
	for s in set(itertools.combinations(list(range(20)), 2)): # generate all variable subset with length of 3
		fout.write('%s %.15f\n' % (''.join([(alphabet[i]) for i in s]), II(list(s), score)))
	t1 = time.time()
	print 'MI time: %d seconds' % (t1-t0)
	'''

	for s in set(itertools.combinations(list(range(20)), 3)): # generate all variable subset with length of 3
		fout.write('%s %.15f\n' % (''.join([(alphabet[i]) for i in s]), deltaN_bar(list(s), score)))
	t2 = time.time()
	print 'DeltaK(3) time: %d seconds' % (t2-t1)

	for s in set(itertools.combinations(list(range(20)), 4)): # generate all variable subset with length of 3
		fout.write('%s %.15f\n' % (''.join([(alphabet[i]) for i in s]), deltaN_bar(list(s), score)))
	t3 = time.time()
	print 'DeltaK(4) time: %d seconds' % (t3-t2)

	pass
if __name__=="__main__":
	main()
