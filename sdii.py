#import numpy.random as np_random
#from scipy.special import binom
from sets import Set
import numpy as np
import itertools
import math
import csv
import time
import sys


__all__=['sdii']

class sdii(object):

	def __init__(self, data):
		self.data = data
		self.entropy_board = {}
		self.weight = np.ones(self.data.shape[0])
		self.isWeighted = False # set True by setWeight()
		self.meff = self.data.shape[0]
		self.varlist = [i for i in xrange(0, self.data.shape[1])]
		self.target = 'all'
		self.order = 2
		self.totalTask = 0


	# currently for msa weight
	def setWeight(self, w):
		self.weight = w
		self.meff = np.sum(w)
		self.isWeighted = True
		print 'set weight vector: %s' % repr(self.weight.shape)
		print 'set meff: %f' % self.meff



	# for task parameters
	def setVarlist(self, varlist):
		self.varlist = varlist
		print 'set varlist: %s ...' % (repr(self.varlist)[0:100])


	# for task parameters
	def setTarget(self, target):
		self.target = target
		print 'set target variable: %s' % self.target


	# for task parameters
	def setOrder(self, order):
		self.order = order
		print 'set order: %d' % self.order


	# for task parameters
	def setTotalTask(self, tn):
		self.totalTask = tn
		print 'set total task number: %d' % self.totalTask

	# weighted entropy test
	'''
	20160304 test ok.
	~/workspace/pdb/test/test_weight/
	data: 
		[[1 1 1 1 1 1]
		 [1 1 1 1 3 3]
		 [3 3 1 1 1 1]
		 [2 1 2 2 3 3]
		 [4 1 2 3 4 4]]

	alphabet sets for column 0 and 1:
		v = [set([1, 2, 3, 4]), set([1, 3])]

	for joint (1, 1) observation 5 rows in [0,1] columns:
		[ 1.  1.  0.  0.  0.]
	weight: [ 0.33333333  0.5         0.5         1.          1.        ]
		p = (1 * 0.3333 + 1 * 0.5 + 0 * 0.5 + 0 * 1 + 0 * 1)/ 3.33 = 0.25
	...
		H = 1.95272419562 vs un-weighted H = 1.92192809489

	'''

	def w_entropy(self, X):
		'''
		#print X.T 
		#print
		H = 0
		#print [set(x) for x in X] 
		#print
		for classes in itertools.product(*[set(x) for x in X]):
		#	print classes
			v = reduce(np.logical_and, (predictions == c for predictions, c in zip(X, classes)))

		#	print [v[i]*self.weight[i] for i in xrange(0, len(self.weight))]
		#	print 
		##	p = np.mean(v) # should divide effective number, which is the sum of all weights
		#	print p
		#	print 
			#p = sum([v[i]*self.weight[i] for i in xrange(0, len(self.weight))])/self.meff
			p = sum(v*self.weight)/self.meff
		#	p = sum([reduce(np.logical_and, (predictions == c for predictions, c in zip(X, classes)))[i]*self.weight[i] for i in xrange(0, len(self.weight))])/self.meff
		#	print p 
		#	print 
			H += -p * np.log2(p) if p > 0 else 0
		return H
		'''
		return np.sum(-p * np.log2(p) if p > 0 else 0 for p in ((sum(reduce(np.logical_and, (predictions == c for predictions, c in zip(X, classes)))*self.weight))/self.meff for classes in itertools.product(*[set(x) for x in X])))

    # general information entropy
    # X: varible set in list type
	def entropy(self, X):
		return np.sum(-p * np.log2(p) if p > 0 else 0 for p in (np.mean(reduce(np.logical_and, (predictions == c for predictions, c in zip(X, classes)))) for classes in itertools.product(*[set(x) for x in X])))

	# calculate entropy with  (-1)^(n+1) sign 
	# for deltaN_bar entropy hashing
	def signed_entropy(self, X, s):
		key = repr(s)
		if key in self.entropy_board:
			H = self.entropy_board[key]
			#print 'found %s' % key
		else:
			if self.isWeighted == False:
				H = math.pow(-1, len(s)+1) * (np.sum(-p * np.log2(p) if p > 0 else 0 for p in (np.mean(reduce(np.logical_and, (predictions == c for predictions, c in zip(X, classes)))) for classes in itertools.product(*[set(x) for x in X]))))
			else:
				H = math.pow(-1, len(s)+1) * self.w_entropy(X)
			self.entropy_board[key] = H
			#print 'put %s' % key
		return H


	def II(self, varset):
		iiv=0.0
		subsets = list(itertools.chain.from_iterable(itertools.combinations(varset, i) for i in range(len(varset)+1)))
		for s in subsets:
			if len(s) == 0:
				continue
			else:
				iiv+=(-1)*self.signed_entropy(self.data[:,s].T, s)
				# non hashing version
				# iiv+=math.pow(-1, len(s))*entropy(data[:,s].T)
		return -iiv


	def deltaN(self, varset):
		dii = {}
		n = len(varset)

		subsets = list(itertools.chain.from_iterable(itertools.combinations(varset, i) for i in range(len(varset)+1)))
		#print 'deltaN_bar()::varset: %s, subsets: %s' % (repr(varset), repr(subsets))
		# varset = Set([2,4,6])
		# [(), (2,), (4,), (6,), (2, 4), (2, 6), (4, 6), (2, 4, 6)]		
		for index in varset:
			deltaX = 0.0
			#subsets = list(itertools.chain.from_iterable(itertools.combinations(varset, i) for i in range(len(varset)+1)))
			for tau in subsets:
				if index in tau:
					# non hashing version
					# deltaX+=math.pow(-1, len(tau)+1)*entropy(data[:,tau].T)
					deltaX+=self.signed_entropy(self.data[:,tau].T, tau) # always hash
			key = '%s\\%d' % (repr(tau), index)
			dii[key] = deltaX
		return dii


	# dependence calculation
	# for a speicific variable set e.g. (X1, X2, X3)
	def deltaN_bar(self, varset):
		deltaX_bar = 1.0
		n = len(varset)

		subsets = list(itertools.chain.from_iterable(itertools.combinations(varset, i) for i in range(len(varset)+1)))
		#print 'deltaN_bar()::varset: %s, subsets: %s' % (repr(varset), repr(subsets))
		# varset = Set([2,4,6])
		# [(), (2,), (4,), (6,), (2, 4), (2, 6), (4, 6), (2, 4, 6)]		
		for index in varset:
			deltaX = 0.0
			#subsets = list(itertools.chain.from_iterable(itertools.combinations(varset, i) for i in range(len(varset)+1)))
			for tau in subsets:
				if index in tau:
					# non hashing version
					# deltaX+=math.pow(-1, len(tau)+1)*entropy(data[:,tau].T)
					deltaX+=self.signed_entropy(self.data[:,tau].T, tau) # always hash
			deltaX_bar*=deltaX

		return math.pow(-1, n)*deltaX_bar


	def sdii_spectrum(self, varset):
		deltaX_bar = 1.0
		n = len(varset)
		entropy_profile = [] # for all the entropy(s)
		deltaX_profile = [] # for all the DeltaX

		subsets = list(itertools.chain.from_iterable(itertools.combinations(varset, i) for i in range(len(varset)+1)))
		#print 'deltaN_bar()::varset: %s, subsets: %s' % (repr(varset), repr(subsets))
		# varset = Set([2,4,6])
		# [(), (2,), (4,), (6,), (2, 4), (2, 6), (4, 6), (2, 4, 6)]		
		for index in varset:
			deltaX = 0.0
			#subsets = list(itertools.chain.from_iterable(itertools.combinations(varset, i) for i in range(len(varset)+1)))
			for tau in subsets:
				if index in tau:
					# non hashing version
					# deltaX+=math.pow(-1, len(tau)+1)*entropy(data[:,tau].T)
					deltaX+=self.signed_entropy(self.data[:,tau].T, tau) # always hash
			deltaX_profile.append(deltaX)
			deltaX_bar*=deltaX

		deltaX_profile.sort() # order does not matter. composition matters

		# store entropy profile
		for i in xrange(1, len(subsets)):
			key = repr(subsets[i])
			entropy_profile.append(self.entropy_board[key]) # order matters

		sdii = math.pow(-1, n)*deltaX_bar

		entropy_profile_str = ','.join([str(v) for v in entropy_profile])
		deltaX_profile_str = ','.join([str(v) for v in deltaX_profile])

		return '%.15f,%s,%s' % (sdii, deltaX_profile_str, entropy_profile_str)




	# integrate deltaN_bar and Interaction information
	def calc_sdii(self, varset):
		n = len(varset)
		if n == 2:
			return self.II(varset)
		else:
			return self.deltaN_bar(varset)


	# threshold calculation
	# for a speicific variable set e.g. (X1, X2, X3)
	def T_l(self, varset):
		deltaX_list = []
		n = len(varset)
		#print 'T_l()::data.shape: %d' % data.shape[0]

		# two variables case we use mutual information instead
		if n==2:
			mi = self.II(varset)
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
						# non hashing version
						# deltaX+=math.pow(-1, len(tau)+1)*entropy(data[:,tau].T)
						deltaX+=self.signed_entropy(self.data[:,tau].T, tau)
			deltaX_list.append(abs(deltaX))

		return math.sqrt(self.data.shape[0]) * min(deltaX_list)




