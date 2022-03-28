import numpy as np
import commp as cp

from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

# input: distance mat & names
def clustering(args):
    assert(len(args)==4), 'Usage: python proc_dd575.py clustering mat.txt name.txt cutoff outfile'

    matfile = args[0]
    namefile = args[1]
    cutoff = float(args[2])
    outfile = args[3]

    mat = np.loadtxt(matfile, delimiter=' ')
    #print(mat.min(), mat.max())
    names = cp.loadlines(namefile)

    dists = squareform(mat)
    z = linkage(dists, 'complete')
    clusters = fcluster(z, t=cutoff, criterion='distance')

    outstr = '\n'.join(['%d %s' % (clusters[i], names[i]) for i in range(len(clusters))])
    with open(outfile, 'w') as fout:
        fout.write('%s\n' % outstr)
    #print(ret)
    cp._info('save %d %f clusters info to %s' % (len(set(clusters)), cutoff, outfile))

if __name__=='__main__':
    cp.dispatch(__name__)