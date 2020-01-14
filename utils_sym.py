import numpy as np
import itertools
import math
import commp as cp

# (weighed) shannon entropy
# varset: varible set in list type 
def _h(data, varset, w='na'):
    X = data[:, varset].T
    if w!='na':
        meff = np.sum(w)
        H = np.sum(-p * np.log2(p) if p > 0 else 0 for p in ((sum(reduce(np.logical_and, (predictions == c for predictions, c in zip(X, classes)))*w))/meff for classes in itertools.product(*[set(x) for x in X])))
    else: # non-weighted ver.
        H = np.sum(-p * np.log2(p) if p > 0 else 0 for p in (np.mean(reduce(np.logical_and, (predictions == c for predictions, c in zip(X, classes)))) for classes in itertools.product(*[set(x) for x in X])))
    return H
    
# for debug purpose
def _wh(data, varset, w):
    X = data[:, varset].T
    meff = np.sum(w)
    print X.T 
    print
    print 'meff: %.f' % meff
    H = 0
    #print [set(x) for x in X] 
    #print
    prob = {}
    for classes in itertools.product(*[set(x) for x in X]):
        print 'class:' + str(classes)
        v = reduce(np.logical_and, (predictions == c for predictions, c in zip(X, classes)))
        print "v:" 
        print v
        # the original p
        p = np.mean(v) # should divide effective number, which is the sum of all weights
        print 'p0: ' + str(p)
        prob[classes] = p

        # apply weights
        p = sum(v*w)/meff
        print 'p1: ' + str(p)
        prob[classes] = p

        # calculate plogp
        H += -p * np.log2(p) if p > 0 else 0
    print
    print '=================='
    print 'probability:'
    print ' '.join(['%s: %s' % (str(k), str(prob[k])) for k in prob])
    print 'sum: %.2f' % (sum(prob.values()))
    return H
#	return np.sum(-p * np.log2(p) if p > 0 else 0 for p in ((sum(reduce(np.logical_and, (predictions == c for predictions, c in zip(X, classes)))*self.weight))/self.meff for classes in itertools.product(*[set(x) for x in X])))

# interaction information
def _ii(data, varset, w='na'):
    iiv=0.0
    subsets = list(itertools.chain.from_iterable(itertools.combinations(varset, i) for i in range(len(varset)+1)))
    # skip the empty set
    for i in xrange(1, len(subsets)):
        s = subsets[i]
        iiv+=pow(-1, len(subsets[i]))*_h(data,subsets[i], w)
    return -iiv


# given varfile output entropy
def mventropy(arglist):


# visualize hierarchical linkage
# for testing weighed entropy
# python utils_sym.py t_linkage t_entropy.data
def t_linkage(arglist):
    import matplotlib.pyplot as plt
    from scipy.cluster.hierarchy import dendrogram
    from scipy.cluster.hierarchy import linkage
    from scipy.cluster.hierarchy import fcluster
    
    data = np.loadtxt(arglist[0], delimiter=',')
    print data.shape
    print data
    max_d = 0.4 # at which distance cutoff we want to cut, 1 - 70%

    linkage_matrix = linkage(data, "complete", metric='hamming')
    print 'linkage matrix:'
    print linkage_matrix
    print
    clusters = fcluster(linkage_matrix, max_d, criterion='distance')
    # clusters contains all the cluster label in the order of data positions
    print clusters
    normdict = cp.freq(clusters)
    print [(1.0/normdict[k]) for k in clusters] # corresponding to the order of data

    # dendrogram
    plt.figure(figsize=(8, 8))
    xt = ['%d-%s' % (i,str(data[i])) for i in xrange(0, len(data))]
    ddata = dendrogram(linkage_matrix, labels=xt, leaf_rotation=90, color_threshold=1, link_color_func=lambda x: "k")
    for i, d, c in zip(ddata['icoord'], ddata['dcoord'], ddata['color_list']):
        x = 0.5 * sum(i[1:3])
        y = d[1]
        plt.plot(x, y, 'o', c=c)
        plt.annotate("%.3g" % y, (x, y), xytext=(0, -5),
                        textcoords='offset points',
                        va='top', ha='center')
        plt.axhline(y=max_d, c='c')   # ['g', 'r', 'c', 'm', 'y', 'k'] 

    plt.tight_layout() # avoid xlabel cutoff
    plt.show()	        
    #plt.savefig(outfile)

# test weighted entropy
def t_wh(arglist):
    data = np.loadtxt(arglist[0], delimiter=',')
    print data.shape
    print data
    
    w = np.loadtxt(arglist[1], delimiter=',')
    print 'weight:'
    print w
    varset=[0]
    print _h(data, varset, w)

# test routine
def t(arglist):
    # python utils_sym.py t t_entropy.data t_entropy.weight
    #
    data = np.loadtxt(arglist[0], delimiter=',')
    w = np.loadtxt(arglist[1], delimiter=',')
    print data.shape
    print data
    varset=[0]
    print _h(data, varset)
    print _h(data, varset,w)

    varset = [0,2]
    print data[:,varset]
    print _h(data, varset, w)
    print _ii(data, varset)
    print _ii(data, varset, w)

    '''
    output:
    (5, 6)
    [[ 1.  1.  1.  1.  1.  1.]
    [ 1.  1.  1.  1.  3.  3.]
    [ 3.  3.  1.  1.  1.  1.]
    [ 2.  1.  2.  2.  3.  3.]
    [ 4.  1.  2.  3.  4.  4.]]
    1.92192809489

    [[ 1.  1.]
    [ 1.  1.]
    [ 3.  1.]
    [ 2.  2.]
    [ 4.  2.]]

    1.92192809489
    0.970950594455

    '''



    '''
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
    '''


if __name__=='__main__':
    cp.dispatch(__name__)