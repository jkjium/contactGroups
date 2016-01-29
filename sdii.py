#!/usr/bin/python
import numpy.random as np_random
from scipy.special import binom
from sets import Set
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
	#n = len(varset)
	#subsets = list(itertools.chain(*[itertools.combinations(range(n), i) for i in range(n+1)]))
	subsets = list(itertools.chain.from_iterable(itertools.combinations(varset, i) for i in range(len(varset)+1)))
	for s in subsets:
		#print s
		if len(s) == 0:
			continue
		else:
			#print data[:,s].T
			iiv+=math.pow(-1, len(s))*entropy(data[:,s].T)
	return -iiv

# dependence calculation
# for a speicific variable set e.g. (X1, X2, X3)
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


# integrate deltaN_bar and Interaction information
def calc_sdii(varset, data):
	n = len(varset)
	print 'calc_sdii(): # of variable: %d' % n
	if n == 2:
		return II(varset, data)
	else:
		return deltaN_bar(varset, data)


# bootstrap on data
def bootstrap(data):
	n = data.shape[0]
	index = np_random.choice(n, n, replace=True)
	#print index
	#print data
	return data[index,:]


# threshold calculation
# for a speicific variable set e.g. (X1, X2, X3)
def T_l(varset, data):
	deltaX_list = []
	n = len(varset)

	# two variables case we use mutual information instead
	if n==2:
		mi = II(varset, data)
		return math.sqrt(data.shape[0])*mi

	for index in varset:
		deltaX = 0.0
		subsets = list(itertools.chain.from_iterable(itertools.combinations(varset, i) for i in range(len(varset)+1)))
		# varset = Set([2,4,6])
		# [(), (2,), (4,), (6,), (2, 4), (2, 6), (4, 6), (2, 4, 6)]
		#print 'T_l()::subsets:'print subsets
		for tau in subsets:
			if index in tau:
				deltaX+=math.pow(-1, len(tau)+1)*entropy(data[:,tau].T)
		deltaX_list.append(abs(deltaX))

	#print deltaX_list
	return math.sqrt(data.shape[0])*min(deltaX_list)


# indicator function on bootstrap data*
# 	count how many T_l >= user input threshold
# t: user threshold 
def G_Bn(data, t, varset, order):
	B = 5 # number of bootstrap
	count = 0
	for i in range(B):
		data_1 = bootstrap(data) # sample with replacement 
		for s in set(itertools.combinations(varset, order)): # generate all variable subset with length of 2
		# varset = Set([2,4,6]), order = 2
		# set([(2, 6), (2, 4), (4, 6)])
			if T_l(list(s), data_1) >= t:
					count+=1
	#print count
	return (1.0/B)*count


# max(1:pk if Tl>=t count++, 1)
def max_Tl_1(data, t, varset, order):
	count = 0 
	for s in set(itertools.combinations(varset, order)): # generate all variable subset with length of 2
	# varset = Set([2,4,6]), order = 2
	# set([(2, 6), (2, 4), (4, 6)])
		if T_l(list(s), data) >= t:
				count+=1
	if count < 1:
		return 1
	else:
		return count*1.0


# find threshold with boostrap
def threshold_t_B(data, alpha, varset, order):
	sk = 4
	pk = binom(len(alphabet), len(varset))
	top = 2*math.sqrt(sk*math.log(pk))
	#print (len(alphabet), sk, pk, top)

	final_t = 0.0
	min_diff = sys.float_info.max

	# get inf(t<=alpha) from all t 
	for t in np.linspace(0,top,10):
		v_G = G_Bn(data, t, varset, order)
		v_m = max_Tl_1(data, t, varset, order)
		#print (v_G, v_m, v_G/v_m, t)
		#rate = G_Bn(data, t, varnum)/max_Tl_1(data, t, varnum)
		diff = alpha - (v_G/v_m)
		if diff > 0 and diff < min_diff:
			min_diff = diff
			final_t = t

	if final_t == 0:
		final_t = top

	return final_t


# forward selection procedure
# return a set of significant variables (index)
def forward_selection(data, alpha, varset, order):
	global alphabet
	ret_varset = Set()
	outfile = 'sdii_%d.txt' % order
	fout = open(outfile, 'w')
	print 'forward_selection()::varset: %s, order: %d' % (repr(varset), order)

	th = threshold_t_B(data, alpha, varset, order)
	print 'forward_selection()::threshold of order [%d]: %f' % (order, th)
	# generate all variable subset with length of order from varlist
	for s in set(itertools.combinations(varset, order)):
		ss = Set(s)
		#print 'forward_selection()::s: %s' % repr(s)
		if len(ss.intersection(varset)) == 0:
			print 'forward_selection()::%s is NOT in %s. skip' % (repr(ss), repr(varset))
			continue
		
		sdii_value = calc_sdii(list(s), data)
		fout.write('%s %.15f\n' % (''.join([(alphabet[i]) for i in s]), sdii_value))
		if sdii_value >= th:
			for var in s:
				ret_varset.add(var)

	print 'forward_selection()::Writing %s' % outfile
	fout.close()
	return ret_varset
		



#alphabet = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
alphabet = ['A','B','C','D']

def main():
	global alphabet

	if len(sys.argv) < 2:
		print 'Usage: python entropy score_file'
		return

	scorefile = sys.argv[1]
	print 'score file: %s' % scorefile

	score = np.loadtxt(scorefile, delimiter=',')
	#print score.shape[0]

	varset = alphabet
	for i in xrange(2,6):
		varset_next = forward_selection(score, 0.05, varset, i)
		if len(varset_next) == 0:
			print 'Main()::stop main loop at order [%d]' % i
			break


	# test forward_selection
	#print forward_selection(score, 0.05, Set([0, 2, 3]), 2)
	

	# test threshold_t_B
	'''
	print threshold_t_B(score, 0.05, 2)
	'''

	# test max_Tl_1
	'''
	print max_Tl_1(score, 5, 3)
	'''
	
	# test T_l
	'''
	for s in set(itertools.combinations(list(range(4)), 3)): # generate all variable subset with length of 2
		print '%s %.15f\n' % (''.join([(alphabet[i]) for i in s]), T_l(list(s), score))
	'''

	# test G_Bn
	# set B = 1
	'''
	G_Bn(score, 0, 2)
	'''

	# test bootstrap
	'''
	bootstrap(score)
	'''
	return

	#print 'score:\n', score
	#print 'AD: %f' % II([2,3], score)
	'''
	alphabet = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
	fout = open(scorefile+'.sdii', 'w')
	print 'calculating mutual information ...'
	t0 = time.time()
	for s in set(itertools.combinations(list(range(20)), 2)): # generate all variable subset with length of 2
		fout.write('%s %.15f\n' % (''.join([(alphabet[i]) for i in s]), II(list(s), score)))

	t1 = time.time()
	print 'MI time: %d seconds' % (t1-t0)

	print 'calculating DeltaK(3) ...'
	for s in set(itertools.combinations(list(range(20)), 3)): # generate all variable subset with length of 3
		fout.write('%s %.15f\n' % (''.join([(alphabet[i]) for i in s]), deltaN_bar(list(s), score)))
	t2 = time.time()
	print 'DeltaK(3) time: %d seconds' % (t2-t1)

	print 'calculating DeltaK(4) ...'
	for s in set(itertools.combinations(list(range(20)), 4)): # generate all variable subset with length of 4
		fout.write('%s %.15f\n' % (''.join([(alphabet[i]) for i in s]), deltaN_bar(list(s), score)))
	t3 = time.time()
	print 'DeltaK(4) time: %d seconds' % (t3-t2)
	'''

if __name__=="__main__":
	main()
