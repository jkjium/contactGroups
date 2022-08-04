import commp as cp
import numpy as np
from sdii import sdii
import itertools



def dtc3(args):
    assert len(args) == 2, 'Usage: python proc_sdiitest.py dtc entropy.mat outifle'
    hmat = np.loadtxt(args[0])
    outfile =args[1]
    # conditional entropy 
    # [(1,), (2,), (3,), (1, 2), (1, 3), (2, 3), (1, 2, 3)]
    m_ch1 = np.array([0,0,0,0,0,-1,1])
    m_ch2 = np.array([0,0,0,0,-1,0,1])
    m_ch3 = np.array([0,0,0,-1,0,0,1])
    
    ch1 = (hmat * m_ch1).sum(axis=1)
    ch2 = (hmat * m_ch2).sum(axis=1)
    ch3 = (hmat * m_ch3).sum(axis=1)
    
    dtc3 = hmat[:,6] - (ch1+ch2+ch3)
    
    np.savetxt(outfile, dtc3.T, fmt='%.6f')
    cp._info('save dct3 vec to %s' % outfile)   


def hxmall(args):
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

    # output: [(1, 2), (1, 3), (2, 3), (1, 2, 3), d3, d2, d1,t(1,2,3)]
    # interaction information
    imat_12  = (hmat * m_h2i_12).sum(axis=1)
    imat_13  = (hmat * m_h2i_13).sum(axis=1)
    imat_23  = (hmat * m_h2i_23).sum(axis=1)
    imat_123 = (hmat * m_h2i_123).sum(axis=1)
    # total correlation
    gmat_123 = (hmat * m_h2g_123).sum(axis=1)

    # differential interaction information
    dii_3 = imat_123 - imat_12
    dii_2 = imat_123 - imat_13
    dii_1 = imat_123 - imat_23
    
    # dual_total correlation
    m_ch1 = np.array([0,0,0,0,0,-1,1])
    m_ch2 = np.array([0,0,0,0,-1,0,1])
    m_ch3 = np.array([0,0,0,-1,0,0,1])
    
    ch1 = (hmat * m_ch1).sum(axis=1)
    ch2 = (hmat * m_ch2).sum(axis=1)
    ch3 = (hmat * m_ch3).sum(axis=1)
    
    dtc3 = hmat[:,6] - (ch1+ch2+ch3)

    imat = np.vstack((imat_12,imat_13,imat_23,imat_123, dii_3, dii_2, dii_1, gmat_123, dtc3))
    np.savetxt(outfile, imat.T, fmt='%.6f')
    cp._info('save II to %s' % outfile)

# kjia@kjia-PC ~/workspace/src 2022-05-02 00:05:14
# $ awk '{print $1,$2,$3}' t.hlist > t.hlist.idx
# $ awk '{print $4,$5,$6,$7,$8,$9,$10}' t.hlist > t.hlist.mat
# https://www.mdpi.com/1099-4300/21/1/88/htm


# entropy to interaction information
# entropy to totoal correlation (multi-information)
# output format: MI_12, MI_13, MI_23, MI_123, DII_3, DII_2, DII_1, TC_123 
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

    # output: [(1, 2), (1, 3), (2, 3), (1, 2, 3), d3, d2, d1,t(1,2,3)]
    imat_12  = (hmat * m_h2i_12).sum(axis=1)
    imat_13  = (hmat * m_h2i_13).sum(axis=1)
    imat_23  = (hmat * m_h2i_23).sum(axis=1)
    imat_123 = (hmat * m_h2i_123).sum(axis=1)
    gmat_123 = (hmat * m_h2g_123).sum(axis=1)

    dii_3 = imat_123 - imat_12
    dii_2 = imat_123 - imat_13
    dii_1 = imat_123 - imat_23

    #print(imat_123, imat_12)

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
# python proc_sdiitest.py entropy_mat t.score.rna
def hmat_test(args):
    assert len(args) == 2, 'Usage: python proc_sdiitest.py entropy_mat_test scorefile order'
    datafile = args[0]
    order = int(args[1])
    #outfile = args[1]
    d = np.loadtxt(datafile, delimiter=',')
    ncol = d.shape[1]
    s = sdii(d)
    s.isWeighted = False
    # generate variable set in the order of :
    # 1,2,3,12,13,23,123 \n
    # 1,2,4,12,14,24,124 \n
    # ...
    varset = list(range(ncol))
    vargroup = list(itertools.chain.from_iterable([itertools.combinations(varset, order)]))
    print('%s choose %d:' % (' '.join(['%d' % i for i in varset]), order))
    print(vargroup)
    print
    retlist = []
    for vs in vargroup: # iterate variable triplets
        idx = ' '.join(['%d' % i for i in vs])
        # generate power set
        spectrum_idx_list = list(itertools.chain.from_iterable(itertools.combinations(vs, i) for i in range(1,order+1)))
        # calculate(hash) entropy of each set in the order of 1,2,3,12,13,23,123
        #entropy_spectrum = map(s.hashed_entropy, spectrum_idx_list)
        #retlist.append('%s %s' % (idx, ' '.join(['%.6f' % e for e in entropy_spectrum])))
        # debug print out
        print(idx)
        print(spectrum_idx_list)
        #print(entropy_spectrum)
    for e in retlist:
        print(e)
    # output file
    with open('t.hlist', 'w') as fout:
        fout.write('%s\n' % '\n'.join(retlist))
    cp._info('save to t.hlist')


# calculate average dii 
# input: 100k.hlist.idx, 100k.hlist.multi
# validated by: 
# paste -d " " 100k.hlist.idx 100k.hlist.multi | \
# awk '{total+=($8+$9+$10);c[$3]++;d[$3]+=$8;c[$2]++;d[$2]+=$9;c[$1]++;d[$1]+=$10} END {for(a in d){print a,d[a],c[a],d[a]/c[a],total/NR;break}}'
def averagedii(args):
    assert len(args) == 3, 'Usage: python proc_sdiitest.py averagedii 100k.hlist.idx 100k.hlist.multi outfile'
    from collections import defaultdict
    idxfile = args[0]
    mfile = args[1]
    outfile = args[2]

    # 100k.hlist.idx
    # 1 2 3
    idxlist = cp.loadtuples(idxfile)
    # 100k.hlist.multi
    # MI_12, MI_13, MI_23, MI_123, DII_3, DII_2, DII_1, TC_123
    mscores = cp.loadtuples(mfile)
    n = len(mscores)
    diisumdict = defaultdict(float)
    idxcount = defaultdict(int) # validate dii numbers
    for i in range(n):
        diisumdict[idxlist[i][2]]+=float(mscores[i][4]) #ii123 - mi12
        diisumdict[idxlist[i][1]]+=float(mscores[i][5]) #ii123 - mi13
        diisumdict[idxlist[i][0]]+=float(mscores[i][6]) #ii123 - mi23

        idxcount[idxlist[i][0]]+=1
        idxcount[idxlist[i][1]]+=1
        idxcount[idxlist[i][2]]+=1

    freqs = set(idxcount.values())
    cp._info(repr(freqs))
    if len(freqs)!=1:
        cp._info('Error: index frequency is no unique')
        return
    else:
        freq = freqs.pop()

    # overall average dii
    total_average = sum(diisumdict.values())/n
    cp._info('total: %.6f, n: %d' % (sum(diisumdict.values()), n))
    # output with triplet form {average.dii3,average.dii2,average.dii1,total_average} 
    fout = open(outfile, 'w')
    for i in range(n):
        #outstr = '%s %.6f %.6f %.6f %.6f' % (' '.join(idxlist[i]),diisumdict[idxlist[i][2]]/freq, diisumdict[idxlist[i][1]]/freq, diisumdict[idxlist[i][0]]/freq, total_average)
        outstr = '%.6f %.6f %.6f %.6f' % (diisumdict[idxlist[i][2]]/freq, diisumdict[idxlist[i][1]]/freq, diisumdict[idxlist[i][0]]/freq, total_average)
        fout.write('%s\n' % outstr)
    fout.close()
    cp._info('save average dii {average.dii3,average.dii2,average.dii1,total_average} to %s' % outfile)


# calculate mean dist of triplets
# distfile:
# python utils_protein2.py writeresdists 1a2t_a.pdb 1a2t_a.pdb.dist
# A 1 A 3 5.5
# index file: {id1, id2, id3}
def mdist(args):
    sdiifile = args[0]
    distfile = args[1]
    outfile = args[2]

    distd = dict(('%s %s' % (d[1],d[3]), float(d[4])) for d in cp.loadtuples(distfile))

    fout=open(outfile, 'w')
    for s in cp.loadtuples(sdiifile):
        k1 = '%s %s' % (s[0], s[1]) if int(s[0]) < int(s[1]) else '%s %s' % (s[1], s[0])
        k2 = '%s %s' % (s[1], s[2]) if int(s[1]) < int(s[2]) else '%s %s' % (s[2], s[1])
        k3 = '%s %s' % (s[0], s[2]) if int(s[0]) < int(s[2]) else '%s %s' % (s[2], s[0])
        md = (distd[k1]+distd[k2]+distd[k3])/3.0
        fout.write('%s %s %s %.8f\n' % (s[0], s[1], s[2], md))
    fout.close()
    cp._info('append average id1,id2,id3,dist to %s' % outfile)





# test sdii.II() and sdii.deltaN_bar() {sdii} with II from entropy operators (hxm)
# $ python proc_sdiitest.py test_ii t.score.rna
# sdii.ii correctness confirmed
# sdii.deltaN_bar() SDII correctness confirmed
def test_ii(args):
    d = np.loadtxt(args[0], delimiter=',')
    s = sdii(d)
    s.isWeighted = False
    varset = [0,1,2,3]
    triple_var = list(itertools.chain.from_iterable([itertools.combinations(varset, 3)]))
    for v3 in triple_var: # iterate variable triplets
        #print(v3, s.II(v3))
        print(v3, s.deltaN_bar(v3))


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