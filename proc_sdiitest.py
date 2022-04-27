import commp as cp
import numpy as np
from sdii import sdii
import itertools

# output matrix form of the entropy spectrum
def entropy_mat(args):
    d = np.loadtxt('t.score', delimiter=',')
    s = sdii(d)
    s.isWeighted = False
    # generate variable set in the order of :
    # 1,2,3,12,13,23,123 \n
    # 1,2,4,12,14,24,124 \n
    # ...
    varset = [0,1,2]
    triple_var = list(itertools.chain.from_iterable([itertools.combinations(varset, 3)]))
    print('all triplets: \n', triple_var)
    for v3 in triple_var:
        print(v3)
        spectrum_idx_list = list(itertools.chain.from_iterable(itertools.combinations(v3, i) for i in range(1,4)))
        print(spectrum_idx_list)
        entropy_spectrum = map(s.hashed_entropy, spectrum_idx_list)
        print(entropy_spectrum)
    #subsets = list(itertools.chain.from_iterable(itertools.combinations(varset, i) for i in range(len(varset))))
    #print(subsets)


# test entropy functions
# validated by matlab scripts
def test_entropy(args):
    print(args)
    d = np.loadtxt('t.score', delimiter=',')
    x=d[:,0]
    y=d[:,1]
    z=d[:,2]
    s = sdii(d)
    s.isWeighted = False
    print('x',s.entropy([x]))
    print('y',s.entropy([y]))
    print('z',s.entropy([z]))
    print('xy',s.entropy([x,y]))
    print('xz',s.entropy([x,z]))
    print('xyz',s.entropy([x,y,z]))
    #print(s.II([0,1]))
    #print(s.II([0,1,2]))
    '''
    print('------------- signed entropy -------------')
    print('xs',s.signed_entropy([x], 'x'))
    print('ys',s.signed_entropy([y], 'y'))
    print('zs',s.signed_entropy([z], 'z'))
    print('xys',s.signed_entropy([x,y], 'xy'))
    print('xzs',s.signed_entropy([x,z], 'xz'))
    '''

    # hashed_entropy
    print('------------- hashed entropy -------------')
    '''
    print('xh',s.hashed_entropy([x], 'x'))
    print('yh',s.hashed_entropy([y], 'y'))
    print('zh',s.hashed_entropy([z], 'z'))
    print('xyh',s.hashed_entropy([x,y], 'xy'))
    print('xzh',s.hashed_entropy([x,z], 'xz'))
    print('xyzh',s.hashed_entropy([x,y,z], 'xyz'))
    '''
    print('xh',s.hashed_entropy([0]))
    print('yh',s.hashed_entropy([1]))
    print('zh',s.hashed_entropy([2]))
    print('xyh',s.hashed_entropy([0,1]))
    print('xzh',s.hashed_entropy([0,2]))
    print('xyzh',s.hashed_entropy([0,1,2]))

if __name__=='__main__':
    cp.dispatch(__name__)