import numpy as np
import itertools
import math
import commp as cp

def foo(arglist):
    print 'hello world'

# general information entropy
# varset: varible set in list type 
def _h(data, varset):
    X = data[:, varset].T
    return np.sum(-p * np.log2(p) if p > 0 else 0 for p in (np.mean(reduce(np.logical_and, (predictions == c for predictions, c in zip(X, classes)))) for classes in itertools.product(*[set(x) for x in X])))

# interaction information
def _ii(data, varset):
    iiv=0.0
    subsets = list(itertools.chain.from_iterable(itertools.combinations(varset, i) for i in range(len(varset)+1)))
	#for s in subsets:
    # skip the empty set
    for i in xrange(1, len(subsets)):
        s = subsets[i]
        iiv+=pow(-1, len(subsets[i]))*_h(data[:,subsets[i]].T)
    return -iiv


def t(arglist):
    # (joint) entropy
    data = np.loadtxt(arglist[0], delimiter=',')
    print data.shape
    X = [0]
    h0 = _h(data[:,X].T)
    print 'H(0): %.8f\n' % h0
    X = [3]
    h3 = _h(data[:,X].T)
    print 'H(3): %.8f\n' % h3
    print 'ii(3): %.8f\n' % _ii(data, X)
    X = [0,3]
    h03 = _h(data[:,X].T)
    print 'H(0,3): %.8f\n' % h03
    # mutual information
    print 'H(0)+H(3)-H(0,3): %.8f\n' % (h0+h3-h03)
    #Xs = [[0,3]]
    print 'ii(0,3): %.8f\n' % _ii(data, X)



if __name__=='__main__':
    cp.dispatch(__name__)