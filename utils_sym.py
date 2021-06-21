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
'''
$ python utils_sym.py t_wh t_entropy.data t_entropy.weight
(5, 6)
[[ 1.  1.  1.  1.  1.  1.]
 [ 1.  1.  1.  1.  3.  3.]
 [ 3.  3.  1.  1.  1.  1.]
 [ 2.  1.  2.  2.  3.  3.]
 [ 4.  1.  2.  3.  4.  4.]]
weight:
[ 0.5  0.5  1.   1.   1. ]
[[ 1.]
 [ 1.]
 [ 3.]
 [ 2.]
 [ 4.]]

meff: 4
alphabet:(1.0,)
v:
[ True  True False False False]
p0: 0.4
p1: 0.25
alphabet:(2.0,)
v:
[False False False  True False]
p0: 0.2
p1: 0.25
alphabet:(3.0,)
v:
[False False  True False False]
p0: 0.2
p1: 0.25
alphabet:(4.0,)
v:
[False False False False  True]
p0: 0.2
p1: 0.25
==================
probability:
(2.0,): 0.25 (3.0,): 0.25 (1.0,): 0.25 (4.0,): 0.25
sum: 1.00
2.0
'''
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
        print 'alphabet:' + str(classes)
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


# return h(wh):value, (joint) probabilities: dictionary
# prob a dictionary: key: tuple, value: float
#    {(1.0, 2.0): 0.0, (3.0, 2.0): 0.0, (3.0, 1.0): 0.201, (2.0, 1.0): 0.0, (2.0, 2.0): 0.201, 
#        (4.0, 2.0): 0.201, (4.0, 1.0): 0.0, (1.0, 1.0): 0.402}
def _hnp(data, varset, w='na'):
    X = data[:, varset].T
    h = 0
    prob = {}
    for classes in itertools.product(*[set(x) for x in X]):
        v = reduce(np.logical_and, (predictions == c for predictions, c in zip(X, classes)))
        p = np.mean(v) if w =='na' else sum(v*w)/np.sum(w) # np.sum(w) is the meff
        prob[classes] = p
        h += -p * np.log2(p) if p > 0 else 0 # entropy
    return h, prob

# for testing purpose
# test _hnp function
def t_hnp(args):
    '''
    $ python utils_sym.py t_hnp t_entropy.data t_entropy.weight 0,2
    [[ 1.  1.]
    [ 1.  1.]
    [ 3.  1.]
    [ 2.  2.]
    [ 4.  2.]]
    1.92192809489
    {(1.0, 2.0): 0.0, (3.0, 2.0): 0.0, (3.0, 1.0): 0.201, (2.0, 1.0): 0.0, (2.0, 2.0): 0.201, 
        (4.0, 2.0): 0.201, (4.0, 1.0): 0.0, (1.0, 1.0): 0.402}

    '''
    data = np.loadtxt(args[0], delimiter=',')
    w = np.loadtxt(args[1], delimiter=',') 
    varset = [int(i) for i in args[2].split(',')]
    H, prob = _hnp(data,varset, w='na')
    #print data[:,varset]
    #print H
    #print prob


# specific information
'''
DeWeese, M. R., & Meister, M. (1999). How to measure the information gained from one symbol. 
    Network: Computation in Neural Systems, 10(4), 325-340.
'''
# I(X|y_i) = H(X) - H(X|y_i)
# H(X|y_i) = -SUM_x(p(x|y_i)log2(x|y_i)) = -SUM_x(p(x,y_i)/p(y_i) * log2(p(x,y_i)/p(y_i)))
# X: target_var
# Y: contribute_var
def _si(data, target_var, contribute_var, w='na'):
    # calculate H(X)
    hX, pX = _hnp(data, [target_var], w) # marginal probability of the target variable
    hY, pY = _hnp(data, [contribute_var], w) # marginal probability of the contribute variable
    hXY, pXY = _hnp(data, [target_var, contribute_var], w) # joint probability

    print 'hX : %.4f' % hX
    print 'pX: '
    print repr(pX)
    print '-----------------------------'
    
    print 'hY : %.4f' % hY
    print 'pY: '
    print repr(pY)
    print '-----------------------------'

    print 'hXY : %.4f' % hXY
    print 'pXY: '
    print repr(pXY)
    print '-----------------------------'

    iX_y = {} # specific information for each observation y \in Y
    for y in pY: # for each observation in the contribute variable, I(X|y_i) = H(X) - H(X|y_i)
        print '------------observation------------ ' + repr(y)

        hX_yi = 0
        for x in pX: # for each observation in the target variable, H(X|y_i) = -SUM_x(p(x,y_i)/p(y_i) * log2(p(x,y_i)/p(y_i)))
            pxy = pXY[(x[0],y[0])]
            if pxy != 0:
                hX_yi+= -pxy/pY[y] * np.log2(pxy/pY[y])
                print 'pXY[(%s,%s)]: %.4f' % (x[0], y[0], pXY[(x[0],y[0])])
                print 'pY[%s]: %.4f' % (y, pY[y])
                print 'pXY[(x[0],y[0])]/pY[y] = %.4f' % (pXY[(x[0],y[0])]/pY[y])
                print 'np.log2(pXY[(x[0],y[0])]/pY[y]) = %.4f' % (np.log2(pXY[(x[0],y[0])]/pY[y]))
                print

        print 'H(X|y)= %.4f'% hX_yi

        iX_yi = hX - hX_yi # I(X|y_i) = H(X) - H(X|y_i)
        iX_y[y] = iX_yi

        print 'I(X|y_i)= %.4f' % iX_y[y]

    sum_IXy=0
    for yi in iX_y:
        sum_IXy+= pY[yi]*iX_y[yi]
        print 'pY[%s]*iX_y[%s]: %.4f * %.4f = %.4f' %(repr(yi), repr(yi), pY[yi], iX_y[yi], pY[yi]*iX_y[yi])
    print 'sum_IXy= %.4f' % sum_IXy
    print

    for yi in iX_y:
        sum_IXy+= pY[yi]*iX_y[yi]
        print '%d %.4f' %(yi[0], pY[yi]*iX_y[yi])
    print

    return iX_y

# test specific information
'''
$ python utils_sym.py t_si t_entropy.data t_entropy.weight 3,4
[3, 4]
[[ 1.  1.  1.  1.  1.  1.]
 [ 1.  1.  1.  1.  3.  3.]
 [ 3.  3.  1.  1.  1.  1.]
 [ 2.  1.  2.  2.  3.  3.]
 [ 4.  1.  2.  3.  4.  4.]]
[ 0.5  0.5  1.   1.   1. ]
hX : 1.3710
pX:
{(2.0,): 0.20000000000000001, (3.0,): 0.20000000000000001, (1.0,): 0.59999999999999998}
-----------------------------
hY : 1.5219
pY:
{(3.0,): 0.40000000000000002, (1.0,): 0.40000000000000002, (4.0,): 0.20000000000000001}
-----------------------------
hXY : 1.9219
pXY:
{(1.0, 3.0): 0.20000000000000001, (3.0, 3.0): 0.0, (3.0, 1.0): 0.0, (2.0, 1.0): 0.0, (2.0, 4.0): 0.0, (2.0, 3.0): 0.20000000000000001, (1.0, 4.0): 0.0, (3.0, 4.0): 0.20000000000000001, (1.0, 1.0): 0.40000000000000002}
-----------------------------
------------observation------------ (3.0,)
pXY[(2.0,3.0)]: 0.2000
pY[(3.0,)]: 0.4000
pXY[(x[0],y[0])]/pY[y] = 0.5000
np.log2(pXY[(x[0],y[0])]/pY[y]) = -1.0000

pXY[(1.0,3.0)]: 0.2000
pY[(3.0,)]: 0.4000
pXY[(x[0],y[0])]/pY[y] = 0.5000
np.log2(pXY[(x[0],y[0])]/pY[y]) = -1.0000

H(X|y)= 1.0000
I(X|y_i)= 0.3710
------------observation------------ (1.0,)
pXY[(1.0,1.0)]: 0.4000
pY[(1.0,)]: 0.4000
pXY[(x[0],y[0])]/pY[y] = 1.0000
np.log2(pXY[(x[0],y[0])]/pY[y]) = 0.0000

H(X|y)= 0.0000
I(X|y_i)= 1.3710
------------observation------------ (4.0,)
pXY[(3.0,4.0)]: 0.2000
pY[(4.0,)]: 0.2000
pXY[(x[0],y[0])]/pY[y] = 1.0000
np.log2(pXY[(x[0],y[0])]/pY[y]) = 0.0000

H(X|y)= 0.0000
I(X|y_i)= 1.3710
pY[(3.0,)]*iX_y[(3.0,)]: 0.4000 * 0.3710 = 0.1484
pY[(1.0,)]*iX_y[(1.0,)]: 0.4000 * 1.3710 = 0.5484
pY[(4.0,)]*iX_y[(4.0,)]: 0.2000 * 1.3710 = 0.2742
sum_IXy= 0.9710

3 0.1484
1 0.5484
4 0.2742

I_XY= 0.9710

'''
def t_si(args):
    datafile = args[0]
    weightfile = args[1]
    varpair = [int(i) for i in args[2].split(',')]
    if len(varpair) != 2:
        cp._err('invalid var index: %s', repr(varpair))
    data = np.loadtxt(datafile, delimiter=',')
    w = np.loadtxt(weightfile, delimiter=',')

    print varpair
    print data
    print w

    t = varpair[0]
    c = varpair[1]

    iX_y = _si(data, t, c, w='na')

    I = _ii(data, [t,c], w='na')
    print 'I_XY= %.4f' % (I)


# interaction information
def _ii(data, varset, w='na'):
    iiv=0.0
    subsets = list(itertools.chain.from_iterable(itertools.combinations(varset, i) for i in range(len(varset)+1)))
    # skip the empty set
    for i in xrange(1, len(subsets)):
        s = subsets[i]
        iiv+=pow(-1, len(subsets[i]))*_h(data,subsets[i], w)
    return -iiv

# variation of information
# VI = 2*H(X,Y) - H(X) - H(Y) *
#   = H(X,Y) - I(X;Y) 
#   = H(X) + H(Y) - 2*I(X;Y)
def _vi(data, var1, var2, w='na'):
    # H(X,Y)
    hxy = _h(data, [var1,var2], w)
    hx = _h(data, [var1], w)
    hy = _h(data, [var2], w)
    #mi = _ii(data, [var1,var2], w)

    print('hxy: %.4f' % hxy)
    print('hx: %.4f' % hx)
    print('hy: %.4f' % hy)
    #print('mi: %.4f' % mi)
    #return hxy - mi
    return 2 * hxy - hx - hy


# test variation of information
def t_vi(args):
    data = np.loadtxt(args[0], delimiter=',')
    varlist = [int(i) for i in args[1].split(',')]
    print data.shape
    print data
    w = np.loadtxt(args[2], delimiter=',') if args[2]!='na' else 'na'
    vi = _vi(data, varlist[0], varlist[1], w)
    print('vi: %.4f' % vi) 
    

# given varfile output entropy
#def mventropy(arglist):


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
    print _wh(data, varset, w)

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