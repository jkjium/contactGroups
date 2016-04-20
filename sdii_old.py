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


# calculate entropy with  (-1)^(n+1) sign 
# for deltaN_bar entropy hashing
def signed_entropy(X, s):
	global entropy_board
	key = repr(s)
	if key in entropy_board:
		H = entropy_board[key]
	else:
		H = math.pow(-1, len(s)+1) * (np.sum(-p * np.log2(p) if p > 0 else 0 for p in (np.mean(reduce(np.logical_and, (predictions == c for predictions, c in zip(X, classes)))) for classes in itertools.product(*[set(x) for x in X]))))
		entropy_board[key] = H
	return H


def entropy_single(X):
	probs = [np.mean(X == c) for c in set(X)]
	#print probs
	return np.sum(-p * np.log2(p) for p in probs)

def II(varset, data, hashing):
	iiv=0.0
	#n = len(varset)
	#subsets = list(itertools.chain(*[itertools.combinations(range(n), i) for i in range(n+1)]))
	subsets = list(itertools.chain.from_iterable(itertools.combinations(varset, i) for i in range(len(varset)+1)))
	for s in subsets:
		#print s
		if len(s) == 0:
			continue
		elif hashing == True:
			#iiv+=math.pow(-1, len(s))*entropy(data[:,s].T)
			# change to hashed version
			iiv+=(-1)*signed_entropy(data[:,s].T, s)
		else:
			iiv+=math.pow(-1, len(s))*entropy(data[:,s].T)
	return -iiv

# dependence calculation
# for a speicific variable set e.g. (X1, X2, X3)
def deltaN_bar(varset, data):
	deltaX_bar = 1.0
	n = len(varset)
	#print 'deltaN_bar(): varset: \n', varset 
	#print 'deltaN_bar(): data: \n' , data
	#sym_1 = ''
	for index in varset:
		deltaX = 0.0
		subsets = list(itertools.chain.from_iterable(itertools.combinations(varset, i) for i in range(len(varset)+1)))
		#print 'deltaN_bar()::varset: %s, subsets: %s' % (repr(varset), repr(subsets))
		# varset = Set([2,4,6])
		# [(), (2,), (4,), (6,), (2, 4), (2, 6), (4, 6), (2, 4, 6)]		
		#sym_str=''
		for tau in subsets:
			if index in tau:
				#deltaX+=math.pow(-1, len(tau)+1)*entropy(data[:,tau].T)
				deltaX+=signed_entropy(data[:,tau].T, tau) # always hash
		deltaX_bar*=deltaX
	#print 'deltaN_bar = %d %s' % (math.pow(-1, n), sym_1)
	return math.pow(-1, n)*deltaX_bar


# integrate deltaN_bar and Interaction information
def calc_sdii(varset, data):
	n = len(varset)
	if n == 2:
		return II(varset, data, True)
	else:
		return deltaN_bar(varset, data)


# bootstrap on data
def bootstrap(n):
	index = np_random.choice(n, n, replace=True)
	#print 'bootstrap()::data1: %s' % repr(index[0:10])
	#return data[index,:]
	return index


# threshold calculation
# for a speicific variable set e.g. (X1, X2, X3)
def T_l(varset, data, hashing):
	deltaX_list = []
	n = len(varset)
	#print 'T_l()::data.shape: %d' % data.shape[0]

	# two variables case we use mutual information instead
	if n==2:
		mi = II(varset, data, hashing)
		#return math.sqrt(data.shape[0])*mi
		return mi

	for index in varset:
		deltaX = 0.0
		subsets = list(itertools.chain.from_iterable(itertools.combinations(varset, i) for i in range(len(varset)+1)))
		# varset = Set([2,4,6])
		# [(), (2,), (4,), (6,), (2, 4), (2, 6), (4, 6), (2, 4, 6)]
		#print 'T_l()::subsets:'print subsets
		for tau in subsets:
			if index in tau:
				if hashing == True:
				#deltaX+=math.pow(-1, len(tau)+1)*entropy(data[:,tau].T)
				# change to hashed version
					deltaX+=signed_entropy(data[:,tau].T, tau)
				else:
					deltaX+=math.pow(-1, len(tau)+1)*entropy(data[:,tau].T)
		deltaX_list.append(abs(deltaX))

	#print deltaX_list
	return math.sqrt(data.shape[0])*min(deltaX_list)

# bootstrap data generation
# index_list: each elemenet is a sampling with replacement
def bootstrap_data(data, index_list, nvar):
	data1 = data[index_list[0], 0]
	for i in xrange(1, nvar):
		col = data[index_list[i], i]
		data1 = np.c_[data1, col]
	return data1



# indicator function on bootstrap data*
# 	count how many T_l >= user input threshold
# t: user threshold 
def G_Bn(data, bootstrap_indexSet, t, varset, order):
	B = len(bootstrap_indexSet) # number of bootstrap
	count = 0
	for i in range(B):
		#print 'G_Bn()::bootstrap # %d' % i
		t0 = time.time()
		#print 'G_Bn()::bootstrap index[1:10]: %s' % str(bootstrap_index[i][0:10])
		data_1 = bootstrap_data(data, bootstrap_indexSet[i], len(varset))
		'''
		print 'G_Bn()::data_1 shape: %s' % repr(data_1.shape)
		print 'G_Bn()::data_1 : %s' % repr(data_1)
		print
		print 'G_Bn()::data : %s' % repr(data)
		'''
		#exit()
		#data_1 = data[bootstrap_index[i],:] # sample with replacement 
		for s in set(itertools.combinations(varset, order)): # generate all variable subset with length of 2
		# varset = Set([2,4,6]), order = 2
		# set([(2, 6), (2, 4), (4, 6)])
			if T_l(list(s), data_1, False) >= t: # DO NOT hashing for bootstrap data
					count+=1
		t1 = time.time()
		#print 'G_Bn: finished in %d seconds' % (t1-t0)
	print 'G_Bn():: # of T >= t : %d, t: %f, count*(1/B): %f' % (count, t, (1.0/B)*count)
	#print count
	return (1.0/B)*count


# max(1:pk if Tl>=t count++, 1)
def max_Tl_1(data, t, varset, order):
	count = 0 
	for s in set(itertools.combinations(varset, order)): # generate all variable subset with length of 2
	# varset = Set([2,4,6]), order = 2
	# set([(2, 6), (2, 4), (4, 6)])
		T = T_l(list(s), data, True) # hashing for real data
		#print 'max_Tl_1()::T: %f' % T
		if T >= t:
				count+=1

	if count < 1:
		print 'max_Tl_1():: # of T >= t: %d, change to 1' % count
		count = 1.0
		return count

	print 'max_Tl_1():: # of T >= t: %d' % count
	return count*1.0


# find threshold with boostrap
def threshold_t_B(data, alpha, varset, order):
	sk = 4
	pk = binom(len(alphabet), order)
	#top = 2*math.sqrt(sk*math.log(pk))
	top = 1.0 
	#print 'threshold_t_B()::'
	#print (len(alphabet), sk, pk, top)

	final_t = 0.0
	min_diff = sys.float_info.max

	n = data.shape[0]
	B = 300

	# a list of lists
	bootstrap_indexSet = []
	for b in xrange(0,B):
		single_var_idx_list = []
		# column-wised samepling, for ith variable
		for i in xrange(0, len(varset)):
			single_var_idx_list.append(np_random.choice(n, n, replace=True))
		bootstrap_indexSet.append(single_var_idx_list)
		#print 'threshold_t_B()::the %dth bootstrap: %s' % (b, repr(single_var_idx_list))
	
	#for i in xrange(0, B):
	#	bootstrap_index.append(np_random.choice(n, n, replace=True))
	# get inf(t<=alpha) from all t 
#	for t in np.linspace(0,top,10):
	for t in np.linspace(0.1,top,10):
		v_G = G_Bn(data, bootstrap_indexSet, t, varset, order)
		#print 'threshold_t_B()::data: %s' % repr(data[1:10,:])
		v_m = max_Tl_1(data, t, varset, order)
		print 'threshold_t_B():: G/Max_T ratio: %f\n' % (v_G/v_m)
		diff = alpha - (v_G/v_m)
		if diff > 0 and diff < min_diff:
			min_diff = diff
			final_t = t
		#break # test sampling

	if final_t == 0:
		final_t = top

	print 'threshold_t_B()::final t: %f' % final_t
	return final_t


# forward selection procedure
# return a set of significant variables (index)
def forward_selection(data, alpha, varset, order):
	global alphabet
	ret_varset = Set()
	outfile = 'sdii_test_%d.txt' % order
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
		
		#print 'forward_selection()::data: %s' % repr(data[1:10,:])
		sdii_value = calc_sdii(list(s), data) # hasing for real data
		'''
		print 'forward_selection()::data: %s' % repr(data[999,:])
		print 'forward_selection:()::list(s): %s' % repr(list(s))
		print 'forward_selection:()::sdii: %f' % sdii_value
		exit()
		'''
		fout.write('%s %.15f\n' % (''.join([(alphabet[i]) for i in s]), sdii_value))
		if sdii_value >= th:
			for var in s:
				ret_varset.add(var)

	print 'forward_selection()::Writing %s' % outfile
	fout.close()
	return ret_varset
		



#alphabet = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
alphabet = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T']
#alphabet = ['X','Y','Z','U','V','W']
#alphabet = ['X1','X2','X3','X4','X5']
entropy_board = {}

def main():
	global alphabet

	if len(sys.argv) < 2:
		print 'Usage: python sdii2.py score_file'
		return

	scorefile = sys.argv[1]
	print 'score file: %s' % scorefile

	score = np.loadtxt(scorefile, delimiter=',')
	#print score.shape[0]


	t1 = time.time()
	varset = range(len(alphabet))
	varset_next = forward_selection(score, 0.1, varset, 2)
	t2 = time.time()
	print varset_next
	print 'use %d seconds' % (t2 - t1)
	
	'''
	print score[0:10,[0,1]].T
	print entropy(score[:,[0,1]].T)
	print entropy(score[:,[0]].T)
	for s in set(itertools.combinations(list(range(6)), 3)): # generate all variable subset with length of 2
		print '%s %.15f\n' % (''.join([(alphabet[i]) for i in s]), deltaN_bar(list(s), score))
	'''
	'''
	varset = range(len(alphabet))
	varset_next = forward_selection(score, 0.05, varset, 2)
	print varset_next
	'''	

	'''
	varset = range(len(alphabet))
	for i in xrange(2,6):
		varset_next = forward_selection(score, 0.05, varset, i)
		if len(varset_next) == 0:
			print 'Main()::stop main loop at order [%d]' % i
			break
		else:
			varset = varset_next
	'''


	# test forward_selection
	#print forward_selection(score, 0.05, Set([0, 2, 3]), 2)
	
	# test T_l
	'''
	for s in set(itertools.combinations(list(range(4)), 3)): # generate all variable subset with length of 2
		print '%s %.15f\n' % (''.join([(alphabet[i]) for i in s]), T_l(list(s), score))
	'''

	'''
	#alphabet = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
	fout = open(scorefile+'.sdii_hash', 'w')
	print 'calculating mutual information ...'
	t0 = time.time()
	for s in set(itertools.combinations(list(range(len(alphabet))), 2)): # generate all variable subset with length of 2
		fout.write('%s %.15f\n' % (''.join([(alphabet[i]) for i in s]), sdii(list(s), score)))

	t1 = time.time()
	print 'MI time: %d seconds' % (t1-t0)

	print 'calculating DeltaK(3) ...'
	for s in set(itertools.combinations(list(range(len(alphabet))), 3)): # generate all variable subset with length of 3
		fout.write('%s %.15f\n' % (''.join([(alphabet[i]) for i in s]), sdii(list(s), score)))
	t2 = time.time()
	print 'DeltaK(3) time: %d seconds' % (t2-t1)
	'''

	'''
	print 'calculating DeltaK(4) ...'
	for s in set(itertools.combinations(list(range(20)), 4)): # generate all variable subset with length of 4
		fout.write('%s %.15f\n' % (''.join([(alphabet[i]) for i in s]), deltaN_bar(list(s), score)))
	t3 = time.time()
	print 'DeltaK(4) time: %d seconds' % (t3-t2)
	'''

if __name__=="__main__":
	main()
