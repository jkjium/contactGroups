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

from sdii import sdii


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
def G_Bn(sdii_obj, bootstrap_indexSet, t, varset, order):
	B = len(bootstrap_indexSet) # number of bootstrap
	count = 0
	for i in range(B):
		#print 'G_Bn()::bootstrap # %d' % i
		t0 = time.time()
		#print 'G_Bn()::bootstrap index[1:10]: %s' % str(bootstrap_index[i][0:10])
		data_1 = bootstrap_data(sdii_obj.data, bootstrap_indexSet[i], len(varset))
		sdii_bootstrap = sdii(data_1) # new hashing object for new data
		'''
		print 'G_Bn()::data_1 shape: %s' % repr(data_1.shape)
		print 'G_Bn()::data_1 : %s' % repr(data_1)
		print
		print 'G_Bn()::data : %s' % repr(data)
		exit()
		'''
		for s in set(itertools.combinations(varset, order)): # generate all variable subset with length of 2
		# varset = Set([2,4,6]), order = 2
		# set([(2, 6), (2, 4), (4, 6)])
			if sdii_bootstrap.T_l(list(s)) >= t: # using the hash table in sdii_bootstrap 
					count+=1
		t1 = time.time()

	print 'G_Bn():: # of T >= t : %d, t: %f, count*(1/B): %f' % (count, t, (1.0/B)*count)
	return (1.0/B)*count


# max(1:pk if Tl>=t count++, 1)
def max_Tl_1(sdii_obj, t, varset, order):
	count = 0 
	for s in set(itertools.combinations(varset, order)): # generate all variable subset with length of 2
	# varset = Set([2,4,6]), order = 2
	# set([(2, 6), (2, 4), (4, 6)])
		T = sdii_obj.T_l(list(s)) # hashing for real data
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
def threshold_t_B(sdii_obj, alpha, varset, order, B):
	sk = 4
	pk = binom(len(varset), order)
	#top = 2*math.sqrt(sk*math.log(pk))
	top = 1.0 
	#print 'threshold_t_B()::'
	#print (len(alphabet), sk, pk, top)

	final_t = 0.0
	min_diff = sys.float_info.max

	n = sdii_obj.data.shape[0]

	# a list of lists
	bootstrap_indexSet = []
	for b in xrange(0,B):
		single_var_idx_list = []
		# column-wised samepling, for ith variable
		for i in xrange(0, len(varset)):
			single_var_idx_list.append(np_random.choice(n, n, replace=True))
		bootstrap_indexSet.append(single_var_idx_list)
		#print 'threshold_t_B()::the %dth bootstrap: %s' % (b, repr(single_var_idx_list))
	
	# get inf(t<=alpha) from all t 
	#for t in np.linspace(0.1,top,10):
	for t in np.linspace(0.0, 0.0005, 20):
		v_G = G_Bn(sdii_obj, bootstrap_indexSet, t, varset, order)
		#print 'threshold_t_B()::data: %s' % repr(data[1:10,:])
		v_m = max_Tl_1(sdii_obj, t, varset, order)
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
def forward_selection(data, alpha, varset, order, B):
	global alphabet
	ret_varset = Set()
	#outfile = 'result_proc_sdii_test_%d.txt' % order
	#fout = open(outfile, 'w')
	print 'forward_selection()::varset: %s, order: %d' % (repr(varset), order)

	sdii_core = sdii(data)
	th = threshold_t_B(sdii_core, alpha, varset, order, B)
	print 'forward_selection()::threshold of order [%d]: %f' % (order, th)

	# generate all variable subset with length of order from varlist
	for s in set(itertools.combinations(varset, order)):
		ss = Set(s)
		#print 'forward_selection()::s: %s' % repr(s)
		if len(ss.intersection(varset)) == 0:
			print 'forward_selection()::%s is NOT in %s. skip' % (repr(ss), repr(varset))
			continue
		
		sdii_value = sdii_core.calc_sdii(list(s)) 

		#fout.write('%s %.15f\n' % (''.join([(alphabet[i]) for i in s]), sdii_value))
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

def main():
	global alphabet

	if len(sys.argv) < 2:
		print 'Usage: python proc_sdii.py score_file'
		return

	scorefile = sys.argv[1]
	print 'score file: %s' % scorefile
	#outfile = scorefile+'.forward'

	score = np.loadtxt(scorefile, delimiter=',')
	#print score.shape[0]
	t1 = time.time()
	varset = range(len(alphabet))
	varset_next = forward_selection(score, 0.1, varset, 20, 300)
	t2 = time.time()
	print varset_next
	print 'use %d seconds' % (t2 - t1)

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
	'''

	t1 = time.time()
	varset = range(len(alphabet))
	th2 = forward_selection(score, 0.1, varset, 2, 300)
	th3 = forward_selection(score, 0.1, varset, 3, 300)
	t2 = time.time()
	'''

	return

	#alphabet = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
	sdii_core = sdii(score)
	fout = open(outfile, 'w')
	print 'calculating mutual information ...'
	t0 = time.time()
	for s in set(itertools.combinations(list(range(len(alphabet))), 2)): # generate all variable subset with length of 2
		fout.write('%s %.15f\n' % (''.join([(alphabet[i]) for i in s]), sdii_core.calc_sdii(list(s))))

	t1 = time.time()
	print 'MI time: %d seconds' % (t1-t0)

	print 'calculating DeltaK(3) ...'
	for s in set(itertools.combinations(list(range(len(alphabet))), 3)): # generate all variable subset with length of 3
		fout.write('%s %.15f\n' % (''.join([(alphabet[i]) for i in s]), sdii_core.calc_sdii(list(s))))
	t2 = time.time()
	print 'DeltaK(3) time: %d seconds' % (t2-t1)



if __name__=="__main__":
	main()
