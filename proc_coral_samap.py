import commp as cp
import numpy as np
import pandas as pd
import scipy as sp
from itertools import groupby
from scipy.spatial import distance
#import matplotlib.pyplot as plt

try:
    import scanpy as sc
    from samap.mapping import SAMAP, prepare_SAMap_loadings
    from samap.analysis import (get_mapping_scores, GenePairFinder, sankey_plot)
    from samap.utils import save_samap, load_samap
    from samalg import SAM    
    from bokeh.plotting import show
    import holoviews as hv
except ImportError:
    cp._info('ignore absent libraries')

class dotdict(dict):
    """dot.notation access to dictionary attributes"""
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

_cscheme_seurat_dimplot=['#f8766d', '#7cae00', '#00bfc4', '#c77cff', '#e68613', '#0cb702', '#00b8e7', '#ed68ed', '#cd9600', '#00be67', '#00a9ff', '#ff61cc', '#aba300', '#00c19a', '#8494ff', '#ff68a1']


# extract gene name, taxon ID from curl download strings
# - adig prost annotations human and drosophila
# $ awk '{printf "curl \"https://rest.uniprot.org/uniprotkb/search?query=%s&fields=gene_names,organism_id\"\necho\n",$1}' 03.uniprot_ID.tsv > batch_curl_all.sh
# $ split -l 140000 batch_curl_all.sh sub.dl. --additional-suffix=.sh
# $ cat *.out > append_info.json.txt
def append_gn_tax(args):
    assert len(args)==2, 'Usage: python proc_coral_samap.py append_gn_tax append_info.json.txt append_info.uid_tax_gn.tsv'
    import json
    from tqdm import tqdm
    infile = args[0]
    outfile = args[1]

    outlist = []
    jsonline=''
    fp = open(infile,'r')
    content = fp.read()
    for line in tqdm(content.split('{"results":'), desc='Processing'):
        if line=='': continue
        jsonline = '{"results":'+line

        #print('New json: '+jsonline)
        jo = json.loads(jsonline)
        try:
            uid = jo['results'][0]['primaryAccession']
            tax = jo['results'][0]['organism'].get('taxonId', '_NA_')
            gid = jo['results'][0].get('genes', '_NA_')
            gn = gid[0].get('geneName', '_NA_') if gid!='_NA_' else '_NA_'
            gnv = gn.get('value', '_NA_') if gn!='_NA_' else '_NA_'
        except Exception as e:
            print(e)
            print(jsonline)
        outlist.append('%s\t%s\t%s' % (uid, tax, gnv))
        jsonline=''
    fp.close()
    with open(outfile, 'w') as fout:
        fout.write('%s\n' % '\n'.join(outlist))

    cp._info('save to append_info.uid_tax_gn.tsv')

# apc procedure fucntions ---------------------------------------------------------------
# tb: tuples from blast [(ad1, sp1), ...,]
# tp: tuples from prost [(ad1, sp1), ...,]
def _hits_comparison(tb, tp): 
    # get full set of protein ids
    bids = set([t[0] for t in tb])
    pids = set([t[0] for t in tp])
    #all_ids = list(bids.union(pids))
    all_ids = list(bids.intersection(pids))
    all_ids.sort()

    # get [id: set(), ... ]
    # adig-s0001.g2:                                                                                     
    #    {'Hvul_g19544_1', 'Hvul_g24627_1', 'Hvul_g23353_1', 'Hvul_g15694_1', 'Hvul_1_1', 'Hvul_g10035_1', 'Hv
    #    ul_g15458_1', 'Hvul_g24294_1', 'Hvul_g28163_1', 'Hvul_g30764_1', 'Hvul_g30759_1', 'Hvul_g7895_1', 'Hv
    #    ul_2_1', 'Hvul_g11044_1', 'Hvul_g27920_1', 'Hvul_g28153_1', 'Hvul_g15015_1', 'Hvul_g24626_1', 'Hvul_g
    #    30756_1', 'Hvul_g24628_1', 'Hvul_g10034_1', 'Hvul_g22665_1'}

    dd = dict((x[0], set([i[1] for i in x[1]])) for x in groupby(tb, lambda t: t[0]))
    dp = dict((x[0], set([i[1] for i in x[1]])) for x in groupby(tp, lambda t: t[0]))

    # iterate all ids
    return [(i,cp.jaccard(dd.get(i, set([])), dp.get(i,set([])))) for i in all_ids]
    

# 1. blast map
# 2. prost map
def maps_comparison(args):
    assert len(args)==3, 'Usage: python maps_comparsion map.blast.txt map.prost.txt outfile'
    tb = [(t[0],t[1]) for t in cp.loadtuples(args[0], delimiter='\t')]
    tp = [(t[0],t[1]) for t in cp.loadtuples(args[1], delimiter='\t')]
    with open(args[2], 'w') as fp:
        fp.write('\n'.join(['%s %.8f' % (t[0], t[1]) for t in _hits_comparison(tb, tp)]))
    cp._info('save jaccard distances to %s' % args[2])
    

############################################################################################
# apc functions  ---------------------------------------------------------------
############################################################################################
#$ python proc_coral_samap.py ut_gene_corr _corr.txt
def ut_gene_corr(args):
    assert len(args) == 1, 'Usage: python proc_coral_samap.py ut_gene_corr t.txt'
    out=_proc_gene_corr_apc(args[0])
    print(out)

# input: gene correlation flat file
# output: np.array: g1, g2, corr_value, apc_value
# outflat = pc._proc_gene_corr_apc('_corr.txt')
# import matplotlib.pyplot as plt
# plt.scatter(outflat[:,2].astype(float),outflat[:,3].astype(float)
# plt.show()
# no sign of apc
def _proc_gene_corr_apc(corrflat):
    # load flat data
    flatdata = np.genfromtxt(corrflat, delimiter=' ', dtype='str')
    r = flatdata[:,0]
    c = flatdata[:,1]
    v = flatdata[:,2].astype(float)
    # calculate apc in matrix form
    ur, iur = np.unique(r, return_inverse=True) # ur: unique row label; iur: replace each r element with their unique label index
    uc, iuc = np.unique(c, return_inverse=True) # ur: unique row label; iur: replace each r element with their unique label index
    import scipy as sp
    outmat = sp.sparse.coo_matrix((v, (iur, iuc)), shape=(len(ur), len(uc))).tocsr()    
    apcmat = _apc_rc(outmat.A)
    # output flat format for visualization
    #apcflat = np.array(['%.8f %s %s' % (apcmat[i[0], i[1]],ur[i[0]], uc[i[1]]) for i in zip(iur, iuc)])
    apcflat = np.array([apcmat[i[0], i[1]] for i in zip(iur, iuc)])
    outflat = np.hstack((flatdata, apcflat[:,None]))
    return outflat



# bipartite apc procedure
# apc_rc = ((row_sum/col_n) outer (col_sum/row_n))/(total_sum/(col_n*row_n))
# sm (score matrix): numpy array score matrix; (pandas.to_numpy())
# r row by c column score matrix
# import matplotlib.pyplot as plt;plt.figure(figsize=(8,6));plt.scatter(sm.flatten(), apc.flatten());plt.show()
def _apc_rc(sm):
    r,c = sm.shape
    '''
    row_mean = sm.sum(axis=1) / c # individual row mean
    col_mean = sm.sum(axis=0) / r # individual column mean
    apc = np.outer(row_mean, col_mean) / (sm.sum() / (r*c))
    '''
    return np.outer(sm.sum(axis=1)/c, sm.sum(axis=0)/r) / (sm.sum()/(r*c))

# symmetrical matrix apc
# apc_ij  = mean_i * mean_j / mean_overall
# apc = (col_sum / (row_n-1)) outer (col_sum / (row_n-1)) / (total_sum / ((col_n^2-col_n))
# input: sm (score matrix): numpy array score matrix
# validation:
# awk '{print $5,$6,$7,$13,$7-$13}' PF01245-6TPQ_i-r2m2i2-mizp3-dcazp3-mipzp3-dcapzp3-dist.vec19 > PF01245-6TPQ_i.i2.mi.mip.apc.txt
# python utils_ce.py cflat2ccmat PF01245-6TPQ_i.i2.mi.mip.apc.txt 0,1 2 PF01245
# sm = = np.loadtxt('PF01245.ccmat')
# apc = _apc_sym(sm)
# $ grep "0.240928" PF012*.txt   
# 1 2 0.571829 0.330901 0.240928 
# apc values agree with APC column in PF01245-6TPQ_i.i2.mi.mip.apc.txt 
# import matplotlib.pyplot as plt;plt.figure(figsize=(8,6));plt.scatter(mi, apc);plt.show()
def _apc_sym(sm):
    return np.outer(sm.sum(axis=0)/(sm.shape[0]-1), sm.sum(axis=0)/(sm.shape[0]-1)) / (sm.sum()/(np.square(sm.shape[0])-sm.shape[0]))

############################################################################################
# h5ad assembling functions  ---------------------------------------------------------------
############################################################################################

# generate sam h5ad file from expression counts
def _samh5ad(h5file,outfile):
    sam = SAM()
    sam.load_data(h5file)
    sam.preprocess_data(
        sum_norm="cell_median",
        norm="log",
        thresh_low=0.0,
        thresh_high=0.96,
        min_expression=1,
    )
    sam.run(
        preprocessing="StandardScaler",
        npcs=100,
        weight_PCs=False,
        k=20,
        n_genes=3000,
        weight_mode='rms'
    )
    sam.leiden_clustering(res=3)
    
    if "PCs_SAMap" not in sam.adata.varm.keys():
        prepare_SAMap_loadings(sam)  

    sam.save_anndata(outfile)
    cp._info('save sam h5ad to %s' % outfile)

# slice adata obj by labels of cluster_name in obs 
# all columms in var will be removed
# h5: input anndata object
# target_obs_column:'Jake'
# labels = ['Jake_10', 'Jake_17', 'Jake_20', 'Jake_22', 'Jake_27', 'Jake_46', 'Jake_54']; corresponds to xe_gastrodermis cell type family
# return: new anndata object contains expression counts for selected cells
def _slice_by_cluster(h5, target_obs_column, labels, obs_droplist=None):
    c = h5.obs[target_obs_column]
    filters = c[c.isin(labels)].index 
    sub_h5 = h5[filters]
    obs = sub_h5.obs.drop(columns=obs_droplist) if obs_droplist is not None else sub_h5.obs
    var = sub_ht.var.drop(columns=sub_h5.var.columns) if sub_h5.var.columns!=0 else sub_h5.var
    return sc.AnnData(X=sub_h5.X, obs=obs, var=var) 
  

# slice data by randomly sample n cells
# n = uniformly sample without replacement (no repeat)
def _slice_by_random(h5, n_cell, obs_droplist):
    filters = np.random.sample(list(range(n_cell)))
    sub_h5 = h5[filters]
    return sc.AnnData(X=sub_h5.X, obs=sub_h5.obs.drop(columns=obs_droplist),  var=sub_h5.var.drop(columns=sub_h5.var.columns))


# append clustering assignments to an h5ad file
# nan will be renamed by unclustered_label
def _append_assignments(h5, clusterfile, cluster_colnames=['cell_type'], index_name='index', unclustered_label='unclustered', delimiter='\t'):
    cluster_assignments = pd.read_csv(clusterfile, header=None, sep=delimiter, index_col=0, names=cluster_colnames)
    cluster_assignments.index.name = index_name 
    h5.obs=h5.obs.merge(cluster_assignments, how='left', left_index=True, right_index=True)    
    for c in cluster_colnames:
        h5.obs[c].fillna(unclustered_label, inplace=True)



# convert RDS converted .mtx {rowname,colname}.tsv files into an h5ad object
# clusterfile: .tsv file without header line
# the h5 and cluster info have the same index name
def _mtx2h5ad(mtxfile, cellname_file, genename_file, index_name='index'):
    import scanpy as sc
    # generating h5ad object
    h5 = sc.read(mtxfile).transpose()
    cellnames = pd.read_csv(cellname_file, header=None, sep='\t', index_col=0)
    cellnames.index.name = index_name
    genenames = pd.read_csv(genename_file, header=None, sep='\t', index_col=0)
    genenames.index.name = index_name
    h5.obs = cellnames
    h5.var = genenames
    #h5.write_h5ad(adig.counts.h5ad.h5ad')
    return h5

# generate h5ad file from mtx file for coral project
# python proc_coral_samap.py h5gen_mtx_cluster adig.counts.mtx adig.colnames.txt adig.rownames.txt ad_cluster_jake.tsv jake adig.counts.c_jake.h5ad
# >>> from proc_coral_samap import _mtx2h5ad, _append_assignments
# >>> h = _mtx2h5ad('adig.counts.mtx', 'adig.colnames.txt', 'adig.rownames.txt')
# >>> _append_assignments(h, 'ad_cluster_jake.tsv', ['jake'])
# h5.write_h5ad('adig.counts.c_jake.h5ad'
def h5gen_mtx_cluster(args):
    assert len(args) == 6, 'Usage: python proc_coral_samap h5adgen mtxfile cellname_file genename_file cluster_file cluster_id outfile'
    h5 = _mtx2h5ad(args[0], args[1], args[2])
    _append_assignments(h5, args[3], [args[4]])
    h5.write_h5ad(args[5])
    cp._info('save to %s' % args[5])



############################################################################################
# alignment score functions --------------------------------------------------------------
############################################################################################

# For generate heatmap
# k=len(np.unique(s.sams['ad'].adata.obs['leiden_clusters']))
# $ ls *prost*family*|awk 'BEGIN {printf "python proc_coral_samap.py combine_scorefiles \""} {printf "%s ", $1} END{printf "\" k outfile"}'
# $ ls *blast*family*|awk 'BEGIN {printf "python proc_coral_samap.py combine_scorefiles \""} {printf "%s ", $1} END{printf "\" k outfile"}'
# outfile: tsv file of combined alignment scores with column ordered by cell type (family) labels
# visualize using the following R code
'''
library(pheatmap)
infile="outfile"
h=pheatmap(t(read.table(infile, row.names=1, header=TRUE, sep='\t')),cluster_rows=F)
source("save_pheatmap.R")
save_pheatmap(h, paste(infile,".png",sep=""), width=12,height=8)
'''
# separator indicate the number of Jake's clusters; all data it is 72; gastrodermis cells: 9 (['Jake_10', 'Jake_17', 'Jake_20', 'Jake_22', 'Jake_27', 'Jake_46', 'Jake_54', |  'Jake_64', 'Jake_25'])
def combine_scorefiles(args):
    assert len(args)==3, 'Usage: python proc_coral_samap.py combine_scores_tsv "list of .tsv files" separator{72} outfile'
    tsvfiles = [fn.strip() for fn in args[0].split(' ') if fn!='']
    cp._info('%s files loaded.' % len(tsvfiles))
    s=int(args[1])
    outfile = args[2]

    #scores = pd.concat([pd.read_csv(fn, sep='\t', index_col=0).iloc[:72,72:] for fn in tsvfiles], axis=1)
    scores = pd.concat([pd.read_csv(fn, sep='\t', index_col=0).iloc[:s,s:] for fn in tsvfiles], axis=1)
    c=list(scores.columns)
    c.sort(key=lambda n: n.split('_')[1:]) # sort by call type families
    scores[c].to_csv(outfile, sep='\t', index=True, header=True)
    cp._info('save all scores %s to %s' % (scores[c].shape, outfile))


############################################################################################
# samap general pipeline functions -------------------------------------------------------
############################################################################################

# call samap main procedure
# sn1,sn2: species names
# sf1,sf2: sam object file names .h5ad (with pre-calculated clustering assignments)
# map_p: homolog graph path
# resolutions: leiden clustering resolution for each species
# assignments: = {"bf":"cluster", "mm":"tissue_celltype"}
# _run_samap_with_h5('ad', '00.adig.counts.Jake_xe_gastodermis.h5ad', 'xe', '00.xesp.counts.Jake_xe_gastodermis.h5ad', samobj=False, map_p='maps.blast/', assignments=assignments)
def _run_samap_with_h5(sn1, fn1, sn2, fn2, samobj=False, map_p='maps/', gnnm=None, resolutions=None, assignments=None, res_dk=1.0, gnnm_transform=True, umap=False):
    from samap.mapping import SAMAP
    from samap.analysis import (get_mapping_scores, GenePairFinder, sankey_plot)
    from samap.utils import save_samap, load_samap
    from samalg import SAM
    import pandas as pd

    if samobj:
        cp._info('loading SAM objects.')
        sam1 = SAM()
        sam1.load_data(fn1)
        sam2 = SAM()
        sam2.load_data(fn2)
        sams = {sn1: sam1, sn2: sam2}
    else:
        cp._info('loading from non-SAM h5 files')
        sams = {sn1:fn1, sn2:fn2}

    cp._info('running SAMap %s%s...' % (sn1,sn2))
    # keys = {"bf":"cluster", "mm":"tissue_celltype"}
    #sm = SAMAP(sams, f_maps=map_p, resolutions=resolutions, keys=assignments)
    sm = SAMAP(sams, f_maps=map_p, gnnm=gnnm, resolutions=resolutions, keys=assignments, res_dk=res_dk, gnnm_transform=gnnm_transform)

    # neigh_from_keys = {'bf':True, 'mm':True}
    #neigh_from_keys = dict((k,True) for k in assignments.keys()) if assignments!=None else None
    if assignments!=None:
        neigh_from_keys={}
        for k in assignments:
            neigh_from_keys[k] = False if assignments[k]=='leiden_clusters' else True
            print(k, neigh_from_keys[k])
    else:
        neigh_from_keys=None
    #sm.run(pairwise=True, neigh_from_keys=neigh_from_keys, umap=False)
    sm.run(pairwise=True, neigh_from_keys=neigh_from_keys, umap=umap) # umap is false for debugging
    #save_samap(sm, 'ad_nv.samap.pkl')
    cp._info('SAMap %s%s done.' %(sn1,sn2))
    return sm

# iterate resolution parameters for leiden_clusters
# comparison between ad_['Jake_10', 'Jake_17', 'Jake_20', 'Jake_22', 'Jake_27', 'Jake_46', 'Jake_54'] and xe_gastrodermis
# seq = np.arange(1,20,4)/10
# for s in seq: _iter_res(s)
def _iter_res(c1='leiden_clusters', c2='leiden_clusters', res1=3, res2=3, output='png'):
    cp._info('resolution iteration: ad: %s res %.2f, xe: %s res %.2f' % (c1, res1, c2, res2))

    assignments={'ad': c1, 'xe': c2}
    resolutions={'ad':res1,'xe':res2}

    sm = _run_samap_with_h5('ad', '00.adig.counts.Jake_xe_gastodermis.h5ad', 'xe', '00.xesp.counts.Jake_xe_gastodermis.h5ad', samobj=False, map_p='maps.blast/', resolutions=resolutions, assignments=assignments)

    D,MappingTable = get_mapping_scores(sm,{'ad':c1,'xe':c2},n_top = 0)
    tsvfile = '02.adxe.xe_gastrodermis.leiden.ad_%.2f_xe_%.2f.tsv' % (res1,res2)
    #MappingTable.to_csv(tsvfile, sep='\t', index=True, header=True)
    #cp._info('save to %s' % tsvfile)

    sk = sankey_plot(MappingTable, align_thr=0.00, species_order = ['ad','xe'])
    hv.save(sk, tsvfile+'.'+output, fmt=output)
    cp._info('save sankey to %s.%s' % (tsvfile, output))
    return (tsvfile, MappingTable)


# prost homology fucntions -------------------------------------------------------------
# convert prost distance to similarity score
# input format: sid1\tsid2\tdistance\te-value
# target_column: 0-based column index, used by pd.iloc
# 
# $ python proc_coral_samap.py prost2similarity 01.prost.adnt.tsv 2 02.prost.adnt.scores.tsv
# $ python proc_coral_samap.py prost2similarity 01.prost.ntad.tsv 2 02.prost.ntad.scores.tsv
# use awk to format it as SAMap mapping
# $ awk -F '\t' '{printf "%s\t%s\tp\tp\tp\tp\tp\tp\tp\tp\t0.0\t%s\n",$1,$2,$5}' 02.prost.adnt.scores.tsv > ad_to_nt.prost.txt
# $ awk -F '\t' '{printf "%s\t%s\tp\tp\tp\tp\tp\tp\tp\tp\t0.0\t%s\n",$1,$2,$5}' 02.prost.ntad.scores.tsv > nt_to_ad.prost.txt
def prost2similarity(args):
    assert len(args) == 3, 'Usage: python proc_coral_samap.py prost2similarity prost.out.tsv target_column ad_to_hy.txt'
    infile = args[0]
    tc = int(args[1])
    outfile = args[2]

    df = pd.read_csv(infile, header=None, sep='\t')
    # ignore all Nan columns
    df = df.dropna(axis=1, how='all')
    # rename previous column names
    df = df.rename(columns=dict((n,i) for i,n in enumerate(df.columns)))

    # append score columns
    i = len(df.columns)
    # inversed log
    df[i+1]= 1/(np.log(df.iloc[:,tc]+1)+1)
    # inversed sqrt
    df[i+2]= 1/(np.sqrt(df.iloc[:,tc]+1)+1)

    df.to_csv(outfile, sep='\t', index=False, header=False)
    cp._info('apprend %d: log, %d: sqrt scores to %s' % (i+1, i+2, outfile))

# for annnotating proteins
# save the goterm information from the top hit uniprot entry
# infile: species separated prost gohits tsv file
# $ grep -v "GO:" adig.proteome.prost.sp.hits_go.tsv |grep  -i "Drosophila melanogaster" > adig.proteome.prost.sp.hits.drome.tsv
# adig-s0001.g1.t2        Q9VCB1  PTOV1_DROME     Protein PTOV1 homolog   Drosophila melanogaster 5889.5  6.35e-06
def nmax_prosthits(args):
    assert len(args) == 2, 'Usage: python proc_coral_samap.py proc_nmax_gohits adig.prost_gohits.drome.tsv adig.prost_gohits.drome.max1.tsv'
    infile = args[0]
    outfile = args[1]
    #df = pd.read_csv(infile, sep="\t", names=['cell_id', 'uniprot', 'u_name', 'u_desc', 'organism', 'dist', 'prost_e-value'])
    df = pd.read_csv(infile, sep="\t", names=['gene_id', 'uniprot', 'u_desc', 'organism', 'dist'])
    idx_min_score = df.groupby('gene_id')['dist'].idxmin()
    df.loc[idx_min_score].to_csv(outfile, sep="\t", index=False, header=False)
    cp._info('save to %s' % outfile)


##-------------------------------------------------------------------------------------
# gene_align_*: remove abnormal AA alphabets gaps = ['.','-',' ','*']
#             : reformat fasta headers; unique ID + transcription ID
##-------------------------------------------------------------------------------------

def _seq_stat(proteome_file, ab_checklist):
    from tqdm import tqdm
    clean_flag = False
    cp._info('checking for sequences contain: %s' % (' '.join(['[%s]' % c for c in ab_checklist])))
    ab_seq = []
    i=0
    for s in tqdm(cp.fasta_iter(proteome_file)):
        for c in ab_checklist:
            if c in s[1]:
                ab_seq.append('[ab]: %s [header]: %s [seq]: %s' % (c, s[0], s[1]))
        i+=1
    print('%d / %d ab sequences found' % (len(ab_seq), i))
    if len(ab_seq)!=0:
        abfile = proteome_file+'.ab.list'
        cp._info('save ab list to %s' % abfile)
        with open(abfile, 'w') as fout:
            fout.write('%s\n' % '\n'.join(ab_seq))
        clean_flag = True
    return clean_flag


# for adig
# s[0]: adig-s0001.g1.t2
def _ad_parser(slist, arg='adig.proteome.rename.fas'):
    _fn_h = lambda sarr: '%s.%s' % (sarr[0], sarr[1])
    seqs = [[_fn_h(s[0].split('.')), s[0], s[1], len(s[1])] for s in slist]
    df = pd.DataFrame(seqs, columns=['output_id','old_id', 'seq', 'length'])
    df.to_csv(arg+'.stat.list', columns=['output_id', 'old_id', 'length'], sep='\t', header=False, index=False)
    idx = df.groupby(['output_id'])['length'].idxmax()
    return df.loc[idx]

# for new nematostella dataset
# input: cleaned sequences slist from proteome.fas 
# output: pd type longest sequences 
def _nt_parser(slist, arg='genes.LUT.tsv'):
    mapping_df = pd.read_csv(arg, sep='\t', names=['id', 'nv2', 'desc', 'old'])
    m = mapping_df.set_index('nv2')['desc'].to_dict()
    
    # NV2.1 == Nv2g000001001 (transcript 1),  Nv2g000001002 (transcript 2)
    seqs = [[m['NV2.' + str(int(s[0][4:10]))], s[0], s[1], len(s[1])] for s in slist]
    df = pd.DataFrame(seqs, columns=['output_id','old_id', 'seq', 'length'])
    # find longest one amount sequences with identical transcript id
    idx = df.groupby(['output_id'])['length'].idxmax()
    return df.loc[idx]

# Spongilla_proteome_70AA.fasta
# map from "c100000_g1_i1_m.41809" to "c100000-g1"
def _sl_parser(slist, arg='Spongilla_proteome_70AA.fasta'):
    _fn_h = lambda sarr: '%s-%s' % (sarr[0], sarr[1])
    seqs = [[_fn_h(s[0].split('_')), s[0], s[1], len(s[1])] for s in slist]

    df = pd.DataFrame(seqs, columns=['output_id','old_id', 'seq', 'length'])
    df.to_csv(arg+'.stat.list', columns=['output_id', 'old_id', 'length'], sep='\t', header=False, index=False)
    idx = df.groupby(['output_id'])['length'].idxmax()
    return df.loc[idx]

# sm=load_samap('02.adxe.samap.prost.jake_cell_type_family.pkl') 
# sm.sams['xe'].adata.var: xe_Xesp_000001
# _fn_get_var: Xesp_000001.t1 to Xesp_000001
def _xe_parser(slist, arg='xesp.proteome.rename.fas'):
    _fn_get_var = lambda sarr: sarr[0]
    seqs = [[_fn_get_var(s[0].split('.')), s[0], s[1], len(s[1])] for s in slist]

    df = pd.DataFrame(seqs, columns=['output_id','old_id', 'seq', 'length'])
    df.to_csv(arg+'.stat.list', columns=['output_id', 'old_id', 'length'], sep='\t', header=False, index=False)
    idx = df.groupby(['output_id'])['length'].idxmax()
    return df.loc[idx]

# same with xe
# sm = load_samap('02.adhy.samap.prost.jake_cell_type_family.pkl')
# sm.sams['hy'].adata.var: hy_Hvul_g10012_1
# _fn_get_var: Hvul_g10012_1.t1 to hy_Hvul_g10012_1 
def _hy_parser(slist, arg='hvul.proteome.rename.fas'):
    _fn_get_var = lambda sarr: sarr[0]
    seqs = [[_fn_get_var(s[0].split('.')), s[0], s[1], len(s[1])] for s in slist]

    df = pd.DataFrame(seqs, columns=['output_id','old_id', 'seq', 'length'])
    df.to_csv(arg+'.stat.list', columns=['output_id', 'old_id', 'length'], sep='\t', header=False, index=False)
    idx = df.groupby(['output_id'])['length'].idxmax()
    return df.loc[idx]


# no need to change
def _sp_parser(slist, arg='spis.proteome.fas'):
    seqs = [[s[0], s[0], s[1], len(s[1])] for s in slist]
    df = pd.DataFrame(seqs, columns=['output_id','old_id', 'seq', 'length'])
    df.to_csv(arg+'.stat.list', columns=['output_id', 'old_id', 'length'], sep='\t', header=False, index=False)
    idx = df.groupby(['output_id'])['length'].idxmax()
    return df.loc[idx]

# no need to change
def _nv_parser(slist, arg='nvec.proteome.rename.fas'):
    seqs = [[s[0], s[0], s[1], len(s[1])] for s in slist]
    df = pd.DataFrame(seqs, columns=['output_id','old_id', 'seq', 'length'])
    df.to_csv(arg+'.stat.list', columns=['output_id', 'old_id', 'length'], sep='\t', header=False, index=False)
    idx = df.groupby(['output_id'])['length'].idxmax()
    return df.loc[idx]

# process proteome data
# 0. remove abnormal alphabets, save stat to *.ab.list
# 1. map header to andata.obs names
# 2. keep the longest sequence for duplicated transcripts
def process_proteome(args):
    assert len(args) == 3, 'Usage: python proc_coral_samap.py process_proteome proteome.fas parser_id outfile'

    _parser_list = {'nt': _nt_parser, 'ad': _ad_parser, 'sl': _sl_parser, 'xe': _xe_parser, 'hy': _hy_parser, 'sp': _sp_parser, 'nv': _nv_parser}

    proteome_file = args[0]
    _fn_s_parser=_parser_list[args[1]]
    cp._info('Parser %s selected' % _fn_s_parser.__name__)
    outfile = args[2]

    # get sequence statistics
    ab_checklist = set(cp.illab).union(set(cp.gaps))
    _seq_stat(proteome_file, ab_checklist)

    # clean sequence and find the longest amount transcripts
    slist=[]
    for h, s in cp.fasta_iter(proteome_file):
        ts = s.translate(str.maketrans('','',''.join(ab_checklist)))
        slist.append((h,ts))

    # processing header and output processed proteome
    cp._info('Filtering short transcripts')
    df_filtered = _fn_s_parser(slist)
    print(df_filtered)
    #df_filtered.to_csv(proteome_file+'.filtered.list', columns=['output_id', 'old_id', 'length'], sep='\t', header=False, index=False)
    cp._info('save filtered header mapping to %s' % proteome_file+'.stat.list')

    # output the converted proteome
    outlist = ['>%s\n%s' % (df_filtered.loc[i]['output_id'], df_filtered.loc[i]['seq']) for i in df_filtered.index]
    with open(outfile, 'w') as fout:
        fout.write('%s\n' % ('\n'.join(outlist)))
    cp._info('save cleaned proteome to %s' % outfile)
    







# kjia@DESKTOP-L0MMU09 ~/workspace/library/dataset/new.nematostella
# for new nematostella dataset
# proteome
# https://simrbase.stowers.org/files/pub/nematostella/Nvec/genomes/Nvec200/aligned/tcs_v2/20240221/NV2g.20240221.protein.fa
# 
# gene name mapping	
# proteome header format in NV2g.20240221.protein.fa
# NV2t025884002.1
# NV2	dataset identifier
# t	stands for transcript
# 025884	gene number
# 002	transcript number
# .1	version 1
# 
# gene name format in GSE200198_alldata.genes.csv
# NVE25884 from NV2t025884002.1
# NVE + strip(025884) = NVE25884
def gene_align_nv2(args):
    assert len(args) == 2, 'Usage: python proc_coral_samap.py gene_align_nv2 NV2g.20240221.protein.fa 00.nv2.proteome.fa'
    infile = args[0]
    outfile = args[1]

    outlist = []
    for h,s in cp.fasta_iter(infile):
        ts = s.translate(str.maketrans('','',''.join(cp.gaps))) # remove '.' at the end of some sequences
        # h: NV2t025884002.1
        new_h = 'NVE' + str(int(h[4:10])) + '_' + h[-5:]
        outlist.append('>%s\n%s' % (new_h, ts))
    
    with open(outfile, 'w') as fout:
        fout.write('%s\n' % '\n'.join(outlist))
    cp._info('save to %s' % outfile)
   
# kjia@DESKTOP-L0MMU09 ~/workspace/samap_coral/stage.nv_ad 2023-12-16 17:23:01
# $ awk '{print substr($1,1,5)}' nvec.genes.tsv|sort|uniq -c
#   33973 Nvec_
# change fasta header to align with nvec.genes.txt 
# from:">v1g152167" to: ">Nvec_v1g152167"
def gene_align_nvec(args):
    assert len(args) == 2, 'Usage: python proc_coral_samap.py gene_align_nvec xenSp1.proteins.fa xesp.proteome.rename.fa'
    infile = args[0]
    outfile = args[1]

    outlist = []
    for h,s in cp.fasta_iter(infile):
        ts = s.translate(str.maketrans('','',''.join(cp.gaps))) # remove '.' at the end of some sequences
        # h = '>v1g152167'
        th = 'Nvec_%s' % h
        outlist.append('>%s\n%s' % (th,ts))
    
    with open(outfile, 'w') as fout:
        fout.write('%s\n' % '\n'.join(outlist))
    cp._info('save to %s' % outfile)


# kjia@DESKTOP-L0MMU09 ~/workspace/samap_coral/stage.xe_ad 2023-12-15 15:19:35
# $ awk '{print substr($1,1,4)}' *genes.tsv|sort|uniq -c
#  29015 Xesp
#   1795 orph
# Xesp records have the same length
# $ grep Xesp *genes.tsv|awk '{print length($1)}'|sort |uniq -c
#  29015 12
# kjia@DESKTOP-L0MMU09 ~/workspace/samap_coral/stage.xe_ad 2023-12-15 15:21:08
# $ grep ">" xenSp1.proteins.fa | awk '{print length($1)}'|sort|uniq -c
#  28010 23

# change fasta header to align with adig.rownames.txt 
# from:">Xe_029129-T2 Xe_029129" to: "Xesp_029129.t2"
def gene_align_xesp(args):
    assert len(args) == 2, 'Usage: python proc_coral_samap.py gene_align_xesp xenSp1.proteins.fa xesp.proteome.rename.fa'
    infile = args[0]
    outfile = args[1]

    import re
    outlist = []
    for h,s in cp.fasta_iter(infile):
        ts = s.translate(str.maketrans('','',''.join(cp.gaps))) # remove '.' at the end of some sequences
        # h = '>Xe_029168-T1 Xe_029168'
        sarr = re.split(r'[_\-\s]',h) # ['Xe', '029168', 'T1', 'Xe', '029168']
        #print(sarr)
        th = 'Xesp_%s.%s' % (sarr[1], sarr[2].lower())
        outlist.append('>%s\n%s' % (th,ts))
    
    with open(outfile, 'w') as fout:
        fout.write('%s\n' % '\n'.join(outlist))
    cp._info('save to %s' % outfile)


# change fasta header to align with adig.rownames.txt 
# from:">Sc4wPfr_165.g10029.t1" to: "Hvul_g10029_1.t1"
def gene_align_hvul(args):
    assert len(args) == 2, 'Usage: python proc_coral_samap.py gene_align_hvul hydra2.0_genemodels.aa hvul.proteome.rename.fasta'
    infile = args[0]
    outfile = args[1]

    outlist = []
    for h,s in cp.fasta_iter(infile):
        ts = s.translate(str.maketrans('','',''.join(cp.gaps))) # remove '.' at the end of some sequences
        #print(">%s\n%s\n" % (h,ts))
        sarr0 = h.split('.') # ['>Sc4wPfr_165', 'g10029', 't1']
        th = 'Hvul_%s_1.%s' % (sarr0[1], sarr0[2]) # Hvul_g10029_1.t1
        outlist.append('>%s\n%s' % (th,ts))
    
    with open(outfile, 'w') as fout:
        fout.write('%s\n' % '\n'.join(outlist))
    cp._info('save to %s' % outfile)


# change fasta header to align with adig.rownames.txt 
# from:">adig_s0001.g1.t2 gene=adig_s0001.g1" to: "adig-s0001.g1.t2"
def gene_align_adig(args):
    assert len(args) == 2, 'Usage: python proc_coral_samap.py gene_align_adig adig.pep.fasta output.fastas'
    infile = args[0]
    outfile = args[1]

    outlist = []
    for h,s in cp.fasta_iter(infile):
        ts = s.translate(str.maketrans('','',''.join(cp.gaps))) # remove '.' at the end of some sequences
        #print(">%s\n%s\n" % (h,ts))
        sarr0 = h.split(' ') # ['>adig_s0001.g1.t2' 'gene=adig_s0001.g1']
        sarr1 = sarr0[0].split('_') # ['adig', 's0001.g1.t2']
        th = '%s-%s' % (sarr1[0], sarr1[1]) # adig-s0001.g1.t2
        outlist.append('>%s\n%s' % (th,ts))
    
    with open(outfile, 'w') as fout:
        fout.write('%s\n' % '\n'.join(outlist))
    cp._info('save to %s' % outfile)

# append meta_cell assignments to ../stage2.blast/01.xxxx.counts_pr.cluster_info.h5ad
# change assignment name to have cell type family appended at the end
# slice expression data by labels=['gastrodermis', 'alga-hosting_cells']
# save to h5ad file for samap ppl
def _proc_stage4_prepare_h5(h5file, meta_assignment_file, outfile):
    s=sc.read_h5ad(h5file)
    _append_assignments(s,meta_assignment_file,cluster_colnames=['metacell'])
    # change cluster name
    # s.obs['cell_type']=s.obs['cell_type'].astype(str)+' ('+s.obs['cell_type_family'].astype(str)+')'
    s.obs['metacell']=s.obs['metacell'].astype(str)+' ('+s.obs['cell_type'].astype(str)+')'
    # slice expression data
    target = 'cell_type_family'
    labels=['gastrodermis', 'alga-hosting_cells']
    n=_slice_by_cluster(s, target, labels, ['leiden_clusters'])
    # save to h5 file
    n.write_h5ad(outfile)
    cp._info('save to %s' % outfile)

# append cell type label to metacell labels
def _proc_anno_metacell():
    sam = SAM()
    sam.load_data('01.hvul.sam.gastodermis.all_info.h5ad')
    sam.adata.obs['metacell']=sam.adata.obs['metacell'].astype(str)+' ('+sam.adata.obs['cell_type'].astype(str)+')'
    sam.save_anndata('01.hvul.sam.gastodermis.all_info.h5ad-1')

    sam = SAM()
    sam.load_data('01.nvec.sam.gastodermis.all_info.h5ad')
    sam.adata.obs['metacell']=sam.adata.obs['metacell'].astype(str)+' ('+sam.adata.obs['cell_type'].astype(str)+')'
    sam.save_anndata('01.nvec.sam.gastodermis.all_info.h5ad-1')

    sam = SAM()
    sam.load_data('01.spis.sam.gastodermis.all_info.h5ad')
    sam.adata.obs['metacell']=sam.adata.obs['metacell'].astype(str)+' ('+sam.adata.obs['cell_type'].astype(str)+')'
    sam.save_anndata('01.spis.sam.gastodermis.all_info.h5ad-1')

    sam = SAM()
    sam.load_data('01.xesp.sam.gastodermis.all_info.h5ad')
    sam.adata.obs['metacell']=sam.adata.obs['metacell'].astype(str)+' ('+sam.adata.obs['cell_type'].astype(str)+')'
    sam.save_anndata('01.xesp.sam.gastodermis.all_info.h5ad-1')

# prepend sid to all the clusters 
def _proc_prepend_sid(samap_file):
    sm = load_samap(samap_file)
    for sid in sm.sams.keys():
        for col in sm.sams[sid].adata.obs:
            sm.sams[sid].adata.obs[col]= sid +'_'+ sm.sams[sid].adata.obs[col].astype(str)
    return sm

# append cell_type annotation to leiden clusters
# target_label: 'leiden_clusters'
# anno_label: 'cell_type'
def _proc_append_anno(samap_file, target_label, anno_label):
    # 02.adhy.samap.blast.gastrodermis.leiden.pkl
    sm=load_samap(samap_file)
    pid = samap_file.split('.')[1]
    cp._info('pid: %s' % pid)

    s1 = sm.sams[pid[:2]].adata.obs
    if anno_label in s1:
        s1[target_label]=s1[target_label].astype(str)+' ('+s1[anno_label].astype(str)+')'

    s2 = sm.sams[pid[2:]].adata.obs
    if anno_label in s2:
        s2[target_label]=s2[target_label].astype(str)+' ('+s2[anno_label].astype(str)+')'

    return sm

# alns = dict((pc._proc_calc_alignments(pc._proc_anno_leiden(sf), sf+'.alignment.tsv')) for sf in samap_files)
# alns for plot sankey plot
def _proc_calc_alignments(sm, outfile):
    ks = list(sm.sams.keys())
    ks.sort()
    pid = ''.join(ks)
    assignments = dict((k,'leiden_clusters') for k in sm.sams.keys())
    D,MappingTable = get_mapping_scores(sm,assignments,n_top = -1)
    MappingTable.to_csv(outfile, sep='\t', index=True, header=True)
    return pid, (D, MappingTable)


# for stage3.gastrodermis comparison
# run samap for all ad vs {sp, hy, xe, nv}
# label: {'cell_type', 'metacell', 'cell_type_family'}
# map_opt: {'blast', 'prost'}
def _proc_gastrodermis(map_opt='blast', label='cell_type'):
    maps= 'maps.%s/' % map_opt
    alignments={}

    assignments = {'ad':'Jake', 'sp':label} if label!=None else None
    outprefix = '02.adsp.samap.%s.gastrodermis.jake_%s' % (map_opt, label) if label!=None else '02.adsp.samap.%s.gastrodermis.leiden' % (map_opt)
    sm = _run_samap_with_h5('ad', '01.adig.sam.gastodermis.all_info.h5ad', 'sp', '01.spis.sam.gastodermis.all_info.h5ad', samobj=True, map_p=maps, assignments=assignments)
    save_samap(sm, outprefix+'.pkl')
    if label==None:
        assignments={'ad':'leiden_clusters', 'sp': 'leiden_clusters'}
    D,MappingTable = get_mapping_scores(sm,assignments,n_top = -1)
    MappingTable.to_csv(outprefix+'.alignment_score.tsv', sep='\t', index=True, header=True)
    alignments['adsp'] = (D, MappingTable)

    assignments = {'ad':'Jake', 'hy':label} if label!=None else None
    outprefix = '02.adhy.samap.%s.gastrodermis.jake_%s' % (map_opt, label) if label!=None else '02.adhy.samap.%s.gastrodermis.leiden' % (map_opt)
    sm = _run_samap_with_h5('ad', '01.adig.sam.gastodermis.all_info.h5ad', 'hy', '01.hvul.sam.gastodermis.all_info.h5ad', samobj=True, map_p=maps, assignments=assignments)
    save_samap(sm, outprefix+'.pkl')
    if label==None:
        assignments={'ad':'leiden_clusters', 'hy': 'leiden_clusters'}
    D,MappingTable = get_mapping_scores(sm,assignments,n_top = -1)
    MappingTable.to_csv(outprefix+'.alignment_score.tsv', sep='\t', index=True, header=True)
    alignments['adhy'] = (D, MappingTable)

    assignments = {'ad':'Jake', 'xe':label} if label!=None else None
    outprefix = '02.adxe.samap.%s.gastrodermis.jake_%s' % (map_opt, label) if label!=None else '02.adxe.samap.%s.gastrodermis.leiden' % (map_opt)
    sm = _run_samap_with_h5('ad', '01.adig.sam.gastodermis.all_info.h5ad', 'xe', '01.xesp.sam.gastodermis.all_info.h5ad', samobj=True, map_p=maps, assignments=assignments)
    save_samap(sm, outprefix+'.pkl')
    if label==None:
        assignments={'ad':'leiden_clusters', 'xe': 'leiden_clusters'}
    D,MappingTable = get_mapping_scores(sm,assignments,n_top = -1)
    MappingTable.to_csv(outprefix+'.alignment_score.tsv', sep='\t', index=True, header=True)
    alignments['adxe'] = (D, MappingTable)

    assignments = {'ad':'Jake', 'nv':label} if label!=None else None
    outprefix = '02.adnv.samap.%s.gastrodermis.jake_%s' % (map_opt, label) if label!=None else '02.adnv.samap.%s.gastrodermis.leiden' % (map_opt)
    sm = _run_samap_with_h5('ad', '01.adig.sam.gastodermis.all_info.h5ad', 'nv', '01.nvec.sam.gastodermis.all_info.h5ad', samobj=True, map_p=maps, assignments=assignments)
    save_samap(sm, outprefix+'.pkl')
    if label==None:
        assignments={'ad':'leiden_clusters', 'nv': 'leiden_clusters'}
    D,MappingTable = get_mapping_scores(sm,assignments,n_top = -1)
    MappingTable.to_csv(outprefix+'.alignment_score.tsv', sep='\t', index=True, header=True)
    alignments['adnv'] = (D, MappingTable)

    return alignments


# save: samap.adata.obs["leiden_clusters_dk"] = pd.Categorical(pc._leiden(samap.adata.uns['Dk'],res=3))
def _leiden(X, res=1, method="modularity", seed = 0):
    if not sp.sparse.isspmatrix_csr(X):
        X = sp.sparse.csr_matrix(X)

    import igraph as ig
    import leidenalg

    sources, targets = X.nonzero()
    weights = X[sources, targets]
    if isinstance(weights, np.matrix):
        weights = weights.A1
    g = ig.Graph(directed=True)
    g.add_vertices(X.shape[0])
    g.add_edges(list(zip(sources, targets)))
    try:
        g.es["weight"] = weights
    except BaseException:
        pass

    if method == "significance":
        cl = leidenalg.find_partition(g, leidenalg.SignificanceVertexPartition,seed=seed)
    else:
        cl = leidenalg.find_partition(
            g, leidenalg.RBConfigurationVertexPartition, resolution_parameter=res,seed=seed
        )

    return np.array(cl.membership)

# calculate leiden cluster using given res
# compare new calculated leiden clusters with leiden_clusters_dk
# ARI and AMI are used for cluster similarity comparison
def _leiden_comp(sam, res):
    import sklearn.metrics as sk
    sam.leiden_clustering(res=res)
    # np conversion is not necessary for now
    a = np.array(sam.adata.obs['leiden_clusters'])
    b = np.array(sam.adata.obs['leiden_clusters_dk'])
    #return [sk.adjusted_mutual_info_score(a,b), sk.adjusted_rand_score(a,b)]
    #return (res, sk.adjusted_rand_score(a,b))
    return (res, sk.adjusted_mutual_info_score(a,b))


############################################################################################
# visualization related --------------------------------------------------------------
############################################################################################
def _heatmap(m):
    import matplotlib.pyplot as plt
    import seaborn as sns
    sns.set(font_scale=0.5)
    fig, ax = plt.subplots(figsize=(14,11))
    sns.heatmap(m, xticklabels=True, yticklabels=True, cmap='YlGnBu', linewidths=0.4, linecolor='lightgray')
    plt.show()

# generate cmap for csleiden sankey plot
# two sids have identical cluster number
# df = sm.sams['ad'].adata.obs['leiden_clusters_dk'] # either one of the sid
# sid = 'ad'
# cmap_id = 'glasbey_hv', 'jet'
def _leiden_dk_sankey_cmap(sids, df, cmap_id='jet', color_scheme=None):
    # filter out singlet clusters in adata.obs assignment
    c=df.value_counts()
    singlets=list(c[c<=1].index)
    filtered_assignments = df[~df.isin(singlets)]

    # assign color code to each cluster label
    clusters = np.unique(filtered_assignments)
    from holoviews.plotting.util import process_cmap
    color_bar = np.array(process_cmap(cmap_id, ncolors=len(clusters), categorical=True)) if color_scheme==None else color_scheme

    # generate cmap for samap.sankey()
    # match ids in MappingTable
    cmap={}
    for sid in sids:
        cmap.update(dict((sid+'_'+str(clusters[i]), color_bar[i % len(color_bar)]) for i in range(len(clusters))))
    return cmap	

# df:cluster_assignments, pandas dataframe, sm.sams['ad'].adata.obs['leiden_clusters']
# n=1 by defualt, for singlet clusters
def _filter_small_clusters(df, n=1):
    c=df.value_counts()
    singlets=list(c[c<=n].index)
    return df[~df.isin(singlets)]

# df_focus: cmap will be generated based on the focused clusters
# return cmap: {'xe_cnidocyte': '#7cae00',..., 'ad_17': '#cd9600'}
# M: alignment scores returned from SAMap.analysis::get_mapping_scores()
def _sankey_cmap(sid1, df1_focus, sid2, df2, M=None, cmap_id='glasbey_hv', color_scheme=None):
    from holoviews.plotting.util import process_cmap
     # filter out singlet clusters in adata.obs assignment
    clusters = np.unique(_filter_small_clusters(df1_focus))
    clusters2 = np.unique(_filter_small_clusters(df2))

    color_bar = np.array(process_cmap(cmap_id, ncolors=len(clusters), categorical=True)) if color_scheme==None else color_scheme
    # generate cmap according to the focused sid;  prepend sid to match ids in the MappingTable
    cmap=dict((sid1+'_'+str(clusters[i]), color_bar[i % len(color_bar)]) for i in range(len(clusters)))
    if M is not None: # find max scores from focused clusters to assign the same color
        sid_map = M.idxmax(axis=1).to_dict()
        cmap.update(dict((sid2+'_'+str(clusters2[i]), cmap[sid_map[sid2+'_'+str(clusters2[i])]]) for i in range(len(clusters2))))
    else: # assign non related color
        cmap.update(dict((sid2+'_'+str(clusters2[i]), color_bar[i % len(color_bar)]) for i in range(len(clusters2))))
    return cmap	   

# t-SNE umap plot with color controlled by:
# 1. cmap_id: matplotlib built-in color schemes like 'glasbey_hv', 'jet'; will be ignored when color_scheme is set
# 2. color_scheme: a list of color names and used cyclically
# 3. sankey_cmap: {'xe_cnidocyte': '#7cae00',..., 'ad_17': '#cd9600'}, cmap_id and color_scheme will be ignored
#    use the color scheme calculated based on MappingTable scores (sankey input)
#
# df = sm.sams['ad'].adata.obs['leiden_clusters_dk']
# xumap = sm.sams['ad'].adata.obsm['X_umap']
def _colored_scatter(sid, df, xumap, cmap_id='jet', color_scheme=None, sankey_cmap=None):
    import matplotlib.pyplot as plt
    from holoviews.plotting.util import process_cmap
    c=df.value_counts()
    singlets=list(c[c<=1].index)
    filtered_assignments = df[~df.isin(singlets)]

    # assign color code to each cluster label
    clusters = np.unique(filtered_assignments)
    if sankey_cmap==None:
        color_bar = np.array(process_cmap(cmap_id, ncolors=len(clusters), categorical=True)) if color_scheme==None else color_scheme
        cmap = dict((clusters[i], color_bar[i % len(color_bar)]) for i in range(len(clusters))) 
    else:
        cmap = dict((clusters[i], sankey_cmap[sid+'_'+str(clusters[i])]) for i in range(len(clusters)))

    plt.figure()
    axes = plt.gca()
    dt=xumap[np.array(~df.isin(singlets))]
    axes.scatter(dt[:,0],dt[:,1], s=10, linewidth=0.0, edgecolor='k', c=list(filtered_assignments.map(cmap)))
    return axes
    #plt.show()

# clusters: pandas dataframe of cluster assignments, for example:
# s.adata.obs
# cn_source: assignment name1, 'cell_type'
# cn_target: assignment name2, 'leiden_clusters, if color_match=True then target color is primary color for color matching
# thr: threshold value
# color_match: match source color to target color with max mapping score
def _sankey_cluster_comp(obs, cn_source,cn_target,thr=0.0, color_match=True, color_scheme=_cscheme_seurat_dimplot, **params):
    import matplotlib.pyplot as plt
    # calcualte index overlaps between clusters
    g1=obs.groupby(obs[cn_source]).groups
    g2=obs.groupby(obs[cn_target]).groups
    # cn1_cluster_labels, cn2_cluster_labels, overlap_number
    #overlap_table = np.array([(k1,k2,len(g2[k2].intersection(g1[k1]))) for k1 in g1 for k2 in g2 if len(g2[k2].intersection(g1[k1]))!=0], dtype=[('source','U11'),('target','U11'),('value','<i8')])
    overlap_table = np.array([[k1,k2,len(g2[k2].intersection(g1[k1]))] for k1 in g1 for k2 in g2 if len(g2[k2].intersection(g1[k1]))!=0])
    R=pd.DataFrame(data=overlap_table[:,0:2], columns=['source','target'])
    R['value']=overlap_table[:,2].astype(int)

    # generate cmap
    c1 = np.unique(R['source'])
    c2 = np.unique(R['target'])
    cmap=dict((c2[i], color_scheme[i % len(color_scheme)]) for i in range(len(c2)))
    if color_match==True: # calculate color matching
        full_score_matrix = R.pivot(index='source', columns='target', values='value').fillna(0)
        sid_map = full_score_matrix.idxmax(axis=1).to_dict() # sid_map[source]=target
        # match source color to target
        cmap.update(dict((c1[i], cmap[sid_map[c1[i]]]) for i in range(len(c1))))
    else: # assign non related color
        cmap.update(dict((c1[i], color_bar[i % len(color_scheme)]) for i in range(len(c1))))

    try:
        from holoviews import dim
        #from bokeh.models import Label
        import holoviews as hv
        hv.extension('bokeh',logo=False)
        hv.output(size=100)        
    except:
        raise ImportError('Please install holoviews-samap with `!pip install holoviews-samap`.')

    depth_map=None
    def f(plot,element):
        plot.handles['plot'].sizing_mode='scale_width'    
        plot.handles['plot'].x_range.start = -600    
        plot.handles['plot'].x_range.end = 1500    

    sankey1 = hv.Sankey(R, kdims=["source", "target"])#, vdims=["Value"])

    #cmap = params.get('cmap','Colorblind')
    label_position = params.get('label_position','outer')
    edge_line_width = params.get('edge_line_width',0)
    show_values = params.get('show_values',False)
    node_padding = params.get('node_padding',4)
    node_alpha = params.get('node_alpha',1.0)
    node_width = params.get('node_width',40)
    node_sort = params.get('node_sort',True)
    frame_height = params.get('frame_height',1000)
    frame_width = params.get('frame_width',800)
    bgcolor = params.get('bgcolor','snow')
    apply_ranges = params.get('apply_ranges',True)

    sankey1.opts(cmap=cmap,label_position=label_position, edge_line_width=edge_line_width, show_values=show_values,
                 node_padding=node_padding,node_cmap=depth_map, node_alpha=node_alpha, node_width=node_width,
                 node_sort=node_sort,frame_height=frame_height,frame_width=frame_width,bgcolor=bgcolor,
                 apply_ranges=apply_ranges,hooks=[f])

    return sankey1    

# calculate element overlaps between clusters of two clustering assignments
def _cluster_overlaps(obs, n1,n2):
    g1=obs.groupby(obs[n1]).groups
    g2=obs.groupby(obs[n2]).groups
    # n1_cluster_label, n2_cluster_label, overlap_number
    return np.array([['%s (%s)' % (str(i), n1), '%s (%s)' % (str(j), n2), len(g2[j].intersection(g1[i]))] for i in g1 for j in g2 if len(g2[j].intersection(g1[i]))!=0])


# for multiple comparision
def _sankey_clusters(obs, cluster_names, **params):
    import matplotlib.pyplot as plt
    pairs = [(cluster_names[i],cluster_names[i+1]) for i in range(len(cluster_names)-1)]
    overlap_table=np.vstack([_cluster_overlaps(obs,n1,n2) for n1,n2 in pairs])
    R=pd.DataFrame(data=overlap_table[:,0:2], columns=['source','target'])
    R['value']=overlap_table[:,2].astype(int)    

    try:
        from holoviews import dim
        #from bokeh.models import Label
        import holoviews as hv
        hv.extension('bokeh',logo=False)
        hv.output(size=100)        
    except:
        raise ImportError('Please install holoviews with `!pip install holoviews`.')

    depth_map=None
    def f(plot,element):
        plot.handles['plot'].sizing_mode='scale_width'    
        plot.handles['plot'].x_range.start = -600    
        plot.handles['plot'].x_range.end = 1500    

    sankey1 = hv.Sankey(R, kdims=["source", "target"])#, vdims=["Value"])

    cmap = params.get('cmap','Colorblind')
    label_position = params.get('label_position','outer')
    edge_line_width = params.get('edge_line_width',0)
    show_values = params.get('show_values',False)
    node_padding = params.get('node_padding',4)
    node_alpha = params.get('node_alpha',1.0)
    node_width = params.get('node_width',40)
    node_sort = params.get('node_sort',True)
    frame_height = params.get('frame_height',1200)
    frame_width = params.get('frame_width', 960)
    bgcolor = params.get('bgcolor','snow')
    apply_ranges = params.get('apply_ranges',True)

    sankey1.opts(cmap=cmap,label_position=label_position, edge_line_width=edge_line_width, show_values=show_values,
                 node_padding=node_padding,node_cmap=depth_map, node_alpha=node_alpha, node_width=node_width,
                 node_sort=node_sort,frame_height=frame_height,frame_width=frame_width,bgcolor=bgcolor,
                 apply_ranges=apply_ranges,hooks=[f])

    return sankey1

def foo(args):
    print(args)


if __name__=='__main__':
    cp.dispatch(__name__)
