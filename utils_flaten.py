import commp as cp
import numpy as np
from scipy.stats import skew

def foo(arglist):
    print arglist

# input: column infile, column index(0-based), outfile
# output: .trans.zscore file with values (+zscore) raw zscore-raw trans zscore-trans
# use cubic root to unskew left screwd distribution
def zscoredist(arglist):
    if len(arglist) < 2:
        cp._err('Usage: python utils_flaten.py zscoredist PF00240.rdca 7')
    infile = arglist[0]
    idx = int(arglist[1]) # column index

    rawlist = []
    for line in cp.loadlines(infile):
        sarr = line.split(' ')
        rawlist.append(float(sarr[idx]))
    oldskew = skew(rawlist)
    #translist = [x ** ( 1. / 3 ) for x in rawlist]
    translist = [cp.root(x,3) for x in rawlist]
    newskew = skew(translist)


    # calculate zscore
    nraw = np.array(rawlist)
    ntrans = np.array(translist)

    m = nraw.mean()
    s = nraw.std()
    zscore_raw = [(x - m)/s for x in nraw]

    m = ntrans.mean()
    s = ntrans.std()
    zscore_trans = [(x - m)/s for x in ntrans]

    # save raw, zscore-raw, trans, zscore-trans
    outfile = '%s.%d.trans.zscore' % (infile, int(idx))
    with open(outfile, 'w') as fout:
        for i in xrange(0, len(rawlist)):
            fout.write('%.6f %.6f %.6f %.6f\n' % (nraw[i], zscore_raw[i], ntrans[i], zscore_trans[i]))
    # output infile, raw_skew, trans_skew
    print '%s %d %.2f %.2f' % (outfile, int(idx), oldskew, newskew)

if __name__=='__main__':
    cp.dispatch(__name__)