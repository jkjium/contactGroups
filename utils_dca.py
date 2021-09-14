import commp as cp
import numpy as np

'''
A python version of mfDCA calculation
The original matlab version was developed by Chris Sander lab.
Marks, Debora S., Thomas A. Hopf, and Chris Sander. "Protein structure prediction from sequence variation." Nature biotechnology 30.11 (2012): 1072-1080.
'''

def dca(args):
    assert len(args) == 4, 'Usage: python utils_dca.py dca scorefile colindexfile weightfile outfile'
    scorefile = args[0]
    colfile = args[1]
    weightfile = args[2]
    outfile = args[3]

    # loading reduced sequences
    score = np.loadtxt(scorefile, delimiter=',', dtype =float).astype(int)
    clist = np.loadtxt(colfile, delimiter=',', dtype=int)
    #print(score.shape)

    w = np.loadtxt(weightfile)
    meff = sum(w)
    cp._info('Loading %s done.' % scorefile)
    cp._info('Loading %s done.' % weightfile)

    # initialize memory
    q = int(np.max(score)+1)
    print(q)
    nrow = score.shape[0]
    ncol = score.shape[1]
    pij_true = np.zeros([ncol,ncol,q,q])
    pi_true = np.zeros([ncol, q])

    # > calculate weighted frequencies
    # single frequency 
    for j in range(nrow):
        for i in range(ncol):
            pi_true[i,score[j,i]] += w[j]
    pi_true = pi_true/meff
    #print(pi_true)
    #print '--------'

    # pair frequency
    for l in range(nrow):
        for i in range(ncol-1): # ncol-1 for i+1 to reach the end
            for j in range(i+1, ncol):
                pij_true[i,j,score[l,i],score[l,j]]+=w[l]
                pij_true[j,i,score[l,j],score[l,i]]=pij_true[i,j,score[l,i],score[l,j]]
    pij_true = pij_true/meff

    # pair frequency diagonal elements
    eyeqxq = np.identity(q)
    for i in range(ncol):
        for aai in range(q):
            for aaj in range(q):
                pij_true[i,i,aai,aaj] = pi_true[i,aai] * eyeqxq[aai, aaj]

    # - [Pij, Pi] = with_pc(Pij_true, Pi_true, pseudocount_weight, alignment_width, q)
    pw = 0.5 # pseudocount_weight set by evfold
    pij = (1-pw)*pij_true + pw/q/q*np.ones([ncol, ncol, q, q])
    pi = (1-pw)*pi_true+pw/q*np.ones([ncol, q])
                    
    # > diagonal
    for i in range(ncol):
        for aai in range(q):
            for aaj in range(q):
                pij[i,i,aai,aaj] = (1-pw) * pij_true[i,i,aai,aaj] + pw/q*eyeqxq[aai,aaj]
    #print('-------------------')
    cp._info('Calculating pij,pi done.')

    # - C = Compute_C(Pij, Pi, alignment_width, q);
    mflat = np.zeros([ncol*(q-1), ncol*(q-1)])
    _fid = lambda _i,_aa,_q: (_q-1)*(_i)+_aa 
    #_fid = lambda _i,_aa,_q: (_q-1)*(_i-1)+_aa + 3
    count=0
    for i in range(ncol):
        for j in range(ncol):
            for aai in range(q-1):
                for aaj in range(q-1):
                    mflat[_fid(i,aai,q), _fid(j,aaj,q)] = pij[i,j,aai,aaj] - pi[i,aai]*pi[j,aaj]
                    '''
                    count+=1
                    if count==200:
                        return
                    print(i,j,aai,aaj,_fid(i,aai,q), _fid(j,aaj,q))
                    '''

    # - invC = inv(C);
    #np.savetxt('pvar', mflat, fmt='%.3f')
    invmflat = np.linalg.inv(mflat)
    cp._info('Calculating invC done.')
    cp._info('Enter DI main loop ...')

    total = cp.ncr(ncol, 2)
    count = 0
    fout = open(outfile, 'w')
    for i in range(ncol-1):
        for j in range(i+1,ncol):
    #for i in range(2):
    #    for j in range(i+1,2):
            # [MI_true, ~, ~] = calculate_mi(i, j, Pij_true, Pi_true, q);
            M = _calc_mi(i,j,pij_true,pi_true,q)
            mfw = np.ones([q,q])

            # > W_m = ReturnW(inC, i, j, q);
            rowlist = [_fid(i,iter,q) for iter in range(0,q-1)]
            collist = [_fid(j,iter,q) for iter in range(0,q-1)]
            mfw[0:q-1, 0:q-1] = np.exp(-invmflat[np.ix_(rowlist, collist)]) # 0:q = 0,1,..,q-1
            #print '---------mfw----------'
            #print('mfw')
            #out = -invmflat[np.ix_(rowlist, collist)]
            #np.savetxt('pvar', out, fmt='%.3f')
            di_mf_pc = _calc_di_mu(i, j, mfw, pi, q)
            fout.write('%d %d %d %d %.6f %.6f\n' % (clist[i],clist[j], i, j, M, di_mf_pc))
            count+=1
            if(count%1000==0):
                cp._info('%d/%d pairs calculated ... ' % (count, total))
    fout.close()
    cp._info('save to %s' % outfile)

# > compute_mu() and compute_di
def _calc_di_mu(i,j,W,p,q):
    epsilon = 1e-4
    diff = 1.0
    mu1 = mu2 = np.ones([1,q])/q
    pi=p[i,:]
    pj=p[j,:]
    _normalize = lambda x: x/np.sum(x)
    # calculate mu1 mu2
    while diff > epsilon:
        scra1 = np.dot(mu2, W.T)
        scra2 = np.dot(mu1, W)

        new1 = _normalize(pi/scra1)
        new2 = _normalize(pj/scra2)

        diff = np.max(np.maximum(np.abs(new1-mu1), np.abs(new2-mu2)))
        mu1 = new1
        mu2 = new2

    # compute_di
    tiny = 1.0e-100;
    #pdir = W * (np.dot(mu1.T,mu2))
    pdir = _normalize(W * (np.dot(mu1.T,mu2)))
    pfac = np.outer(p[i,:], p[j,:])

    di = np.trace(np.dot(pdir.T,np.log((pdir+tiny)/(pfac+tiny))))
    return di

# calculate mutual information
def _calc_mi(i, j, p2, p1, q):
    M = 0.0
    for a in range(q):
        for b in range(q):
            if p2[i,j,a,b]>0:
                M+=p2[i,j,a,b]*np.log(p2[i,j,a,b]/p1[i,a]/p1[j,b])
    return M

if __name__ == '__main__':
        cp.dispatch(__name__)