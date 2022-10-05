import commp as cp
import numpy as np

# prepare data for plotting accumulative_accuracy_ratio vs ranking_ratio(top X percentage ranked)
# extended from utils_ps3.py appendrankp
# input: ce value file(.vec15), ce value indeces {1,2,3},  response_index, response_cutoff
# output: .png
# $ python proc_apc3.py rankp RF00001.hlist_3j7o_7.pdb.dist_r3.md.tc3.sdii3.ii3.dtc3.vec15 6 3 5.0 500
def rankp(args):
    assert len(args) == 5, 'Usage: python utils_ce.py rankp cefile collist{1,2,3,} yi yv 500'
    from scipy.stats import rankdata

    infile = args[0]
    cols = [int(c) for c in args[1].split(',')]
    yi = int(args[2])
    yv = float(args[3])
    num_sample = int(args[4])

    #.vec15: r3.md.tc3.sdii3.ii3.dtc3
    cp._info('loading data ...')
    data = np.loadtxt(infile)
    n = len(data) # number of pairs/triplets
    outlist = []
    #for i in cols:
    for i in range(len(cols)):
        # ranking of ith cevalue
        p = rankdata(-data[:,cols[i]], method='min') # get rank value for ith ce column

        # generating 0,1 response then paste with rank values
        # r: 01 rank, .vec15[:,3] is mdist (response y)
        r =np.c_[(data[:,yi]<yv).astype(float), p]

        # sort vec2 by rank values
        # s: vec2: 01, sorted_rank
        s = r[r[:,1].argsort()]

        total_contact = sum(r[:,0]) # sum number of TRUE values
        
        # sL: 01, rank ratio (high to low order)
        #sL = np.c_[s[:,0], s[:,1]/L]
        sL = np.c_[s[:,0], s[:,1]/n] # sL[:,1] in [1/n, 1]

        # calculate accumulative number of contact
        cp._info('calculating accumulative contact numbers')
        a = np.cumsum(sL[:,0])

        # outplt: x-axis: sorted_rank_ratio y-axis:accumulative_accuracy_ratio
        outplt = np.c_[sL[:,1], a/total_contact]

        # sample 500 points and save to ce separated files 
        sample_idx = np.round(np.linspace(0, len(outplt) - 1, num_sample)).astype(int)
        outdata = outplt[sample_idx,1]
        outfile = '%s.%d.ppv' % (infile,cols[i])
        np.savetxt(outfile, outdata, fmt='%.6f', newline = " ")
        cp._info('%s :: save %dth ppv to %s' % (infile, cols[i],outfile))

        '''
        # plot ppv vs sorted_rank for debugging
        cp._info('plotting rank/L vs acc')
        import matplotlib.pyplot as plt
        #plt.plot(outdata[:,0], outdata[:,1], label=legends[i])
        plt.plot(outdata[:,0], outdata[:,1])

    # finishing plots
    plt.xlabel('ranking ratio')
    plt.ylabel('ppv')
    plt.legend(loc='best')
    plt.tight_layout()
    plt.show()
    '''

# $ python proc_apc3.py cerankl RF00001.hlist_3j7o_7.pdb.dist_r3.md.tc3.sdii3.ii3.dtc3.vec15 6 120 out
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