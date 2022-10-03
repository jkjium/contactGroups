import commp as cp
import numpy as np

# prepare data for plotting accumulative_accuracy_ratio vs ranking_ratio(top X percentage ranked)
# extended from utils_ps3.py appendrankp
# input: ce value file(.vec15), ce value indeces {1,2,3},  response_index, response_cutoff
# output: .png
def rankp(args):
    assert len(args) == 6, 'Usage: python utils_ce.py cerankl cefile collist{1,2,3,} yi yv legends outfile'
    from scipy.stats import rankdata
    import matplotlib.pyplot as plt

    infile = args[0]
    cols = [int(c) for c in args[1].split(',')]
    yi = int(args[2])
    yv = float(args[3])
    legends = args[4]
    outfile = args[5]

    #.vec15: r3.md.tc3.sdii3.ii3.dtc3
    cp._info('loading data ...')
    data = np.loadtxt(infile)
    n = len(data) # number of pairs/triplets
    outlist = []
    for i in cols:
        # get rank of ith ce
        cp._info('ranking of %dth cevalue' % i)
        p = rankdata(-data[:,i], method='min') # get rank value for ith ce column

        # generating 0,1 response then paste with rank values
        # r: 01 rank, .vec15[:,3] is mdist (response y)
        r =np.c_[(data[:,yi]<yv).astype(float), p]

        # sort vec2 by rank values
        # s: vec2: 01, sorted_rank
        cp._info('sorting by ranks ...')
        s = r[r[:,1].argsort()]

        total_contact = sum(r[:,0]) # sum of 1s
        cp._info('total contact: %d' % total_contact)
        
        # sL: 01, rank ratio (high to low order)
        #sL = np.c_[s[:,0], s[:,1]/L]
        sL = np.c_[s[:,0], s[:,1]/n] # sL[:,1] in [1/n, 1]
        # masked sorted rank/L
        # msL: 01 sorted_rank/L
        #msL = sL[sL[:,1] <=10.0, :]
        msL = sL

        # calculate accumulative number of contact
        cp._info('calculating accumulative contact numbers')
        a = np.cumsum(msL[:,0])

        # outplt: x-axis: sorted_rank_ratio y-axis:accumulative_accuracy_ratio
        outplt = np.c_[msL[:,1], a/total_contact]

        # save to ce separate files 

        # ppv, sorted_rank
        cp._info('plotting rank/L vs acc')
        plt.plot(outplt[:,0], outplt[:,1])

    # finishing plots
    plt.xlabel('ranking ratio')
    plt.ylabel('ppv')
    plt.legend(loc='best')
    plt.savefig(outfile)    
    plt.tight_layout()
    plt.show()

if __name__=='__main__':
    cp.dispatch(__name__)