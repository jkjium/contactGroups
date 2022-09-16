import commp as cp
import numpy as np

# rank ce value as percentage or x * L
# extended from utils_ps3.py appendrankp
# input: .vec15
# output: .png
def cerankl(args):
    assert len(args) == 4, 'Usage: python utils_ce.py cerankl cefile collist{1,2,3,} L outfile'
    from scipy.stats import rankdata
    import matplotlib.pyplot as plt

    infile = args[0]
    cols = [int(c) for c in args[1].split(',')]
    L = float(args[2])
    outfile = args[3]

    #.vec15: r3.md.tc3.sdii3.ii3.dtc3
    cp._info('loading .vec15 data ...')
    data = np.loadtxt(infile)
    n = len(data)
    outlist = []
    for i in cols:
        # get rank of ith ce
        cp._info('ranking of %dth cevalue' % i)
        p = rankdata(-data[:,i], method='min')
        # r: 01 rank, .vec15[:,3] is mdist
        r =np.c_[(data[:,3]<5).astype(float), p]
        # sort by rank
        # s: 01 sorted_rank
        cp._info('sorting ranks ...')
        s = r[r[:,1].argsort()]
        total_contact = sum(r[:,0])
        cp._info('total contact: %d' % total_contact)
        
        # sL: 01, sorted_rank/L
        #sL = np.c_[s[:,0], s[:,1]/L]
        sL = np.c_[s[:,0], s[:,1]/n]
        # masked sorted rank/L
        # msL: 01 sorted_rank/L
        #msL = sL[sL[:,1] <=10.0, :]
        msL = sL

        # calculate accumulative number of contact
        cp._info('calculating accumulative contact numbers')
        a = np.cumsum(msL[:,0])
        #a = [sum(s[0:j+1,0]) for j in range(len(msL))]

        # outplt: sorted_rank/L acc
        outplt = np.c_[msL[:,1], a/total_contact]
        # ppv = acc[:,2] / acc[:,1]

        # ppv, sorted_rank
        #outplt = np.c_[ppv, acc[:,1]/L]
        #print(outplt)
        cp._info('plotting rank/L vs acc')
        plt.plot(outplt[:,0], outplt[:,1])
        plt.tight_layout()
        plt.show()

    '''
        
    outnplist = np.array(outlist)
    np.savetxt(outfile, outnplist.T, fmt='%.6f')
    cp._info('save L/x rank to %s' % outfile)
    '''

if __name__=='__main__':
    cp.dispatch(__name__)