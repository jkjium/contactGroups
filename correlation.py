import numpy as np
import commp as cp
import itertools
import math

# general information entropy
# X: varible set in list type
def entropy(X):
	return np.sum(-p * np.log2(p) if p > 0 else 0 for p in (np.mean(reduce(np.logical_and, (predictions == c for predictions, c in zip(X, classes)))) for classes in itertools.product(*[set(x) for x in X])))

def w_entropy(X):
	H = 0
	print repr(X)
	print [set(x) for x in X]
	#print [set(x) for x in X] 
	#print
	#for classes in itertools.product(*[set(x) for x in X]):
	#	print classes
	#	v = reduce(np.logical_and, (predictions == c for predictions, c in zip(X, classes)))
	#	print [v[i]*self.weight[i] for i in xrange(0, len(self.weight))]
	#	print 
	##	p = np.mean(v) # should divide effective number, which is the sum of all weights
	#	print p
	#	print 
		#p = sum([v[i]*self.weight[i] for i in xrange(0, len(self.weight))])/self.meff
	#	p = sum(v*self.weight)/self.meff
	#	p = sum([reduce(np.logical_and, (predictions == c for predictions, c in zip(X, classes)))[i]*self.weight[i] for i in xrange(0, len(self.weight))])/self.meff
	#	print p 
	#	print 
	#	H += -p * np.log2(p) if p > 0 else 0
	#return H
	#return np.sum(-p * np.log2(p) if p > 0 else 0 for p in ((sum(reduce(np.logical_and, (predictions == c for predictions, c in zip(X, classes)))*self.weight))/self.meff for classes in itertools.product(*[set(x) for x in X])))

# for testing purpose
def main():
	pass


if __name__ == '__main__':
	main()