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

# print total number of genes in a flat cluster file
def printtotal(args):
    assert len(args) == 1, 'Usage: python proc_dd575.py totalgenes flatclusterfile'
    flatfile = args[0]

    glist = []
    for c in cp.loadtuples(flatfile):
        for g in c[7].split(','):
            glist.append(g)
    cp._info('total genes in %s : %d / %d' % (flatfile, len(glist), len(set(glist))))
    #print ('%s\n' % '\n'.join(glist))

# output list of genes in a flat cluster file    
def flat2vec(args):
    assert len(args) == 2, 'Usage: python proc_dd575.py flat2vec clusterfile outfile'
    flatfile = args[0]
    outfile =args[1]

    glist =  []
    for c in cp.loadtuples(flatfile):
        for g in c[7].split(','):
            glist.append(g)
    with open(outfile, 'w') as fout:
        fout.write('%s\n' % '\n'.join(glist))
    cp._info('save %d genes to %s' % (len(glist), outfile))


# after cutoff expansion
# get rest single genes that cannot form clusters
def findsingle(args):
    assert len(args) == 3, 'Usage: python proc_dd575.py findsingle 12100.f.out total.clusters out.single.list'
    reffile =args[0]
    clusterfile = args[1]
    outfile = args[2]

    rlist = set() 
    for r in cp.loadtuples(reffile): # cluster information using 12100 as cutoff
        for g in r[7].split(','):
            rlist.add(g)
    cp._info('%d genes found in reference clusters' % len(rlist))

    clist = set()
    for c in cp.loadtuples(clusterfile):
        for g in c[7].split(','):
            clist.add(g)
    cp._info('%d genes found in selected clusters' % len(clist))

    outlist = rlist - clist
    with open(outfile, 'w') as fout:
        fout.write('%s\n' %('\n'.join(outlist)))
    cp._info('save %d singlets to %s' % (len(outlist),outfile))


# output cluster information {id, number of members, tm.mean, tm.sd, tm.min, tm.max, member_list}
# input: the output from function clustering: cluster_id, member_id, sorted by cluster_id from func clustering()
# input: pairwise dist/similarity file, in this case tmalign scores 
def flatcluster(args):
    assert len(args) == 4, 'Usage: python proc_dd575.py flatcluster clustering.out tmscore.list filterlist outfile'

    clusterfile = args[0]
    scorefile = args[1]
    filterfile = args[2] # selected clusters with high tmscores
    outfile = args[3]

    # load dist/similarity file into a dictionary
    # AG6000091 AG6000243 0.36034
    # AG6000091 AG6000244 0.40890
    #print(cp.loadtuples(dfile))
    sdict = dict(('%s %s' % (s[0], s[1]), float(s[2])) for s in cp.loadtuples(scorefile))
    #  print(take(10, dinfodict.iteritems())) # python 3.6

    # load clustering information
    # 431 AG6000091
    # 521 AG6000243
    full_cinfo = cp.loadtuples(clusterfile)

    # load filered clusters
    # 15500 7 2 0.6328 0.6328 0.6328 0.0000 AG6000799,AG6032177
    filterlist = []
    for s in cp.loadtuples(filterfile):
        for g in s[7].split(','):
            filterlist.append(g)

    # filter out selected genes
    cinfo = [c for c in full_cinfo if c[1] not in filterlist]

    cp._info('full: %d - select: %d : current: %d' % (len(full_cinfo), len(filterlist), len(cinfo)))

    # ('224', ['AG6006935', 'AG6007576', 'AG6010406', 'AG6033388']) ('xxx', [xxx])
    cdict = dict((k, [i[1] for i in list(g)]) for k, g in itertools.groupby(cinfo, lambda x: x[0]))

    # iterate all clusters to map pairwise scores (tm)
    # for each cluster ID
    fout = open(outfile, 'w')
    for c in cdict:
        genes = cdict[c]
        #print(c, cdict[c])
        # single member cluster
        if len(genes) < 2:
            outstr = '%s %d 0.0 0.0 0.0 0.0 %s' % (c, len(genes), ','.join(genes))
            #print(outstr)
            fout.write('%s\n' % outstr)
            continue
        # multi-member cluster
        scores = []
        for i in range(len(genes)):
            for j in range(i+1, len(genes)):
                k1 = '%s %s' % (genes[i], genes[j])
                k2 = '%s %s' % (genes[j], genes[i])
                score = sdict[k1] if k1 in sdict else sdict[k2]
                scores.append(score)
        ns = np.array(scores)
        #print(ns)
        # cluster_id, cluster_len, tm.mean, tm.min, tm.max, tm.std
        outstr = '%s %d %.4f %.4f %.4f %.4f %s' % (c, len(genes), ns.mean(), ns.min(), ns.max(), ns.std(), ','.join(genes))
        fout.write('%s\n' % outstr)
        #print(outstr)

        # output pymol scripts
        pml = []
        pml.append('delete all')
        for n in genes:
            pml.append('load %s.pdb' % n)
        for n in genes[1:]:
            pml.append('cealign %s, %s\ncenter' % (genes[0], n))
        with open('%s' % c, 'w') as fml:
            fml.write('%s\n' % '\n'.join(pml))
        
    fout.close()
    cp._info('save %d flatclusters to %s' % (len(cdict), outfile))
    cp._info('cluster_id, cluster_len, tm.mean, tm.min, tm.max, tm.std')


# compare two clusters generated using different cutoffs
# input: cluster.12100.out cluster.13000.out (must be sorted by the clusters ID first)
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