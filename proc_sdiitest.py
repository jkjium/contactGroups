import commp as cp
import numpy as np
from sdii import sdii
import itertools

# kjia@kjia-PC ~/workspace/src 2022-05-02 00:05:14
# $ awk '{print $1,$2,$3}' t.hlist > t.hlist.idx
# $ awk '{print $4,$5,$6,$7,$8,$9,$10}' t.hlist > t.hlist.mat
# https://www.mdpi.com/1099-4300/21/1/88/htm


# entropy to interaction information
# entropy to totoal correlation (multi-information)
# output format: MI_12, MI_13, MI_23, MI_123, TC_123 
# 2nd TC is equivalent to MI
def hxm(args):
    assert len(args) == 2, 'Usage: python proc_sdiitest.py hxm entropy.mat outfile'
    hmat = np.loadtxt(args[0])
    outfile = args[1]
    # [(1,), (2,), (3,), (1, 2), (1, 3), (2, 3), (1, 2, 3)]
    m_h2i_123= np.array([1,1,1,-1,-1,-1,1])
    m_h2i_12 = np.array([1,1,0,-1,0,0,0])
    m_h2i_13 = np.array([1,0,1,0,-1,0,0])
    m_h2i_23 = np.array([0,1,1,0,0,-1,0])
    m_h2g_123= np.array([1,1,1,0,0,0,-1])

    #print(hmat)
    #print(hmat * m_h2i_123)

    # output: [(1, 2), (1, 3), (2, 3), (1, 2, 3), t(1,2,3)]
    imat_12  = (hmat * m_h2i_12).sum(axis=1)
    imat_13  = (hmat * m_h2i_13).sum(axis=1)
    imat_23  = (hmat * m_h2i_23).sum(axis=1)
    imat_123 = (hmat * m_h2i_123).sum(axis=1)
    gmat_123 = (hmat * m_h2g_123).sum(axis=1)

    dii_3 = imat_123 - imat_12
    dii_2 = imat_123 - imat_13
    dii_1 = imat_123 - imat_23

    print(imat_123, imat_12)

    imat = np.vstack((imat_12,imat_13,imat_23,imat_123, dii_3, dii_2, dii_1, gmat_123))
    np.savetxt(outfile, imat.T, fmt='%.6f')
    cp._info('save II to %s' % outfile)


# output matrix form of the entropy spectrum
# $ python proc_sdiitest.py entropy_mat
# $ cat t.hlist
# 0 1 2 1.918296 1.459148 0.000000 1.918296 1.918296 1.459148 1.918296
# 0 1 3 1.918296 1.459148 0.650022 1.918296 2.251629 1.792481 2.251629
# 0 2 3 1.918296 0.000000 0.650022 1.918296 2.251629 0.650022 2.251629
# 1 2 3 1.459148 0.000000 0.650022 1.459148 1.792481 0.650022 1.79248
def entropy_mat(args):
    d = np.loadtxt('t.score', delimiter=',')
    s = sdii(d)
    s.isWeighted = False
    # generate variable set in the order of :
    # 1,2,3,12,13,23,123 \n
    # 1,2,4,12,14,24,124 \n
    # ...
    varset = [0,1,2,3]
    triple_var = list(itertools.chain.from_iterable([itertools.combinations(varset, 3)]))
    print('all triplets: \n', triple_var)
    retlist = []
    for v3 in triple_var: # iterate variable triplets
        idx = ' '.join(['%d' % i for i in v3])
        # generate power set
        spectrum_idx_list = list(itertools.chain.from_iterable(itertools.combinations(v3, i) for i in range(1,4)))
        # calculate(hash) entropy of each set in the order of 1,2,3,12,13,23,123
        entropy_spectrum = map(s.hashed_entropy, spectrum_idx_list)
        retlist.append('%s %s' % (idx, ' '.join(['%.6f' % e for e in entropy_spectrum])))
        # debug print out
        print(idx)
        print(spectrum_idx_list)
        print(entropy_spectrum)
    for e in retlist:
        print(e)
    # output file
    with open('t.hlist', 'w') as fout:
        fout.write('%s\n' % '\n'.join(retlist))
    cp._info('save to t.hlist')


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