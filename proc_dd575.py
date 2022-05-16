import numpy as np
import commp as cp
import itertools
import subprocess as sp 

from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

# cmd call tmalign (tmscore.sh)
# input two gene names (in string)
def _tmscore(param):
    return sp.Popen(['tmscore.sh', param[0], param[1]], stdout=sp.PIPE).communicate()[0]

# test _tmscore function
# input two gene names
def tmscore(args):
    assert len(args) == 2, 'Usage: python proc_dd575.py tmscore AG6000091 AG6000243'
    print(_tmscore(args))

# pair all the genes
def pairgenes(args):
    assert len(args) == 2, 'Usage: python proc_dd575.py pairgenes full.list out.pair.list'
    genes = cp.loadlines(args[0])
    fout = open(args[1], 'w')
    for i in range(len(genes)):
        for j in range(i+1, len(genes)):
            fout.write('%s %s\n' % (genes[i], genes[j]))
    cp._info('save paired genes to %s' % args[1])


# compare two clusters generated using different cutoffs
# input: cluster.12100.out cluster.13000.out (must be sorted first)
# cluster.xxx.out format:
# 4 AG6000735
# 23 AG6000741
# 23 AG6000767
def cluster_comp(args):
    assert len(args) ==3, 'Usage: python proc_dd575.py cluster_comp cluster.12100.out cluster.13000.out'
    oldfile = args[0]
    newfile = args[1]
    outfile = args[2]

    c_old = cp.loadtuples(oldfile)
    c_new = cp.loadtuples(newfile)

    # ('224', ['AG6006935', 'AG6007576', 'AG6010406', 'AG6033388']) ('xxx', [xxx])
    c_old_dict = dict((k, [i[1] for i in list(g)]) for k, g in itertools.groupby(c_old, lambda x: x[0]))
    c_new_dict = dict((k, [i[1] for i in list(g)]) for k, g in itertools.groupby(c_new, lambda x: x[0]))

    if len(c_new_dict) > len(c_old_dict):
        cp._err('args order error: old:%s new:%s' % (oldfile, newfile))

    fout = open(outfile, 'w')
    for p in c_old_dict:
        oldc = set(c_old_dict[p])
        for q in c_new_dict:
            newc = set(c_new_dict[q])
            if len(newc)>3:
                c_inter = oldc.intersection(newc)
                diff = newc - oldc
                if len(c_inter)!=0 and len(diff)!=0:
                    if len(oldc)<4:
                        fout.write('oldc: %s n: %d newc: %s n: %d diff: %s\n' % (p, len(oldc), q, len(newc), ','.join([g for g in newc])))
                    else:
                        fout.write('oldc: %s n: %d newc: %s n: %d diff: %s\n' % (p, len(oldc), q, len(newc), ','.join([g for g in diff])))

    fout.close()
    cp._info('save comparison result to %s' % outfile)


# input: distance mat & names
# output sorted cluster information:
# 4 AG6000735
# 23 AG6000741
# 23 AG6000767
# python proc_dd575.py clustering 575.s.pairwise.pkl.txt 575.names.pkl.txt 7200 out
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

    c_sort = sorted([(clusters[i], names[i]) for i in range(len(clusters))], key = lambda x: x[0])
    outstr = '\n'.join(['%d %s' % (c[0], c[1]) for c in c_sort])
    #outstr = '\n'.join(['%d %s' % (clusters[i], names[i]) for i in range(len(clusters))])
    with open(outfile, 'w') as fout:
        fout.write('%s\n' % outstr)
    #print(ret)
    cp._info('save %d %f clusters info to %s' % (len(set(clusters)), cutoff, outfile))

if __name__=='__main__':
    cp.dispatch(__name__)