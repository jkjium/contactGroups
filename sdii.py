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
		else:
			H = math.pow(-1, len(s)+1) * (np.sum(-p * np.log2(p) if p > 0 else 0 for p in (np.mean(reduce(np.logical_and, (predictions == c for predictions, c in zip(X, classes)))) for classes in itertools.product(*[set(x) for x in X]))))
			self.entropy_board[key] = H
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


	# dependence calculation
	# for a speicific variable set e.g. (X1, X2, X3)
	def deltaN_bar(self, varset):
		deltaX_bar = 1.0
		n = len(varset)

		for index in varset:
			deltaX = 0.0
			subsets = list(itertools.chain.from_iterable(itertools.combinations(varset, i) for i in range(len(varset)+1)))
			#print 'deltaN_bar()::varset: %s, subsets: %s' % (repr(varset), repr(subsets))
			# varset = Set([2,4,6])
			# [(), (2,), (4,), (6,), (2, 4), (2, 6), (4, 6), (2, 4, 6)]		
			for tau in subsets:
				if index in tau:
					# non hashing version
					# deltaX+=math.pow(-1, len(tau)+1)*entropy(data[:,tau].T)
					deltaX+=self.signed_entropy(self.data[:,tau].T, tau) # always hash
			deltaX_bar*=deltaX

		return math.pow(-1, n)*deltaX_bar


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



