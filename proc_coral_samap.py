import commp as cp
import numpy as np
import pandas as pd
from itertools import groupby
from scipy.spatial import distance

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
    

# bipartite apc procedure
# apc_rc = ((row_sum/col_n) outer (col_sum/row_n))/(total_sum/(col_n*row_n))
# sm: numpy array score matrix; (pandas.to_numpy())
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

# apc procedure for symmetrical matrix
# apc_ij  = mean_i * mean_j / mean_overall
# apc = (col_sum / (row_n-1)) outer (col_sum / (row_n-1)) / (total_sum / ((col_n^2-col_n))
# input: sm: numpy array score matrix
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
def _slice_by_cluster(h5, target_obs_column, labels, obs_droplist):
    c = h5.obs[target_obs_column]
    filters = c[c.isin(labels)].index 
    sub_h5 = h5[filters]
    return sc.AnnData(X=sub_h5.X, obs=sub_h5.obs.drop(columns=obs_droplist),  var=sub_h5.var.drop(columns=sub_h5.var.columns))


# append clustering assignments to an h5ad file
# nan will be renamed by unclustered_label
def _append_assignments(h5, clusterfile, cluster_colnames=['cell_type'], index_name='index', unclustered_label='unclustered'):
    cluster_assignments = pd.read_csv(clusterfile, header=None, sep='\t', index_col=0, names=cluster_colnames)
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
def _run_samap_with_h5(sn1, fn1, sn2, fn2, samobj=False, map_p='maps/', resolutions=None, assignments=None):
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
    sm = SAMAP(sams, f_maps=map_p, resolutions=resolutions, keys=assignments)

    # neigh_from_keys = {'bf':True, 'mm':True}
    neigh_from_keys = dict((k,True) for k in assignments.keys()) if assignments!=None else None
    #sm.run(pairwise=True, neigh_from_keys=neigh_from_keys, umap=False)
    sm.run(pairwise=True, neigh_from_keys=neigh_from_keys, umap=False) # umap is false for debugging
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


# homology graph fucntions -------------------------------------------------------------
# !! incomplete
# convert prost output to blast fmt6 format
# adig-s0001.g1.t2        Hvul_1_1.g28296                         5523.5  1.12e-10
# ['adig-s0001.g1.t2', 'Hvul_1_1.g28296', '', '', '', '5523.5', '1.12e-10'] 
# blast -outfmt 6 # https://www.metagenomics.wiki/tools/blast/blastn-output-format-6
# 1      2      3      4      5        6       7      8     9     10    11    12
# qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
def prost2fmt6(args):
    assert len(args) == 2, 'Usage: python proc_coral_samap.py prost2fmt6 prost.out.tsv ad_to_hy.txt'
    infile = args[0]
    outfile = args[1]

    max_d = 0.0
    for t in cp.loadtuples(infile, delimiter='\t'):
        d = float(t[5])
        if d > max_d: max_d = d
    print(max_d)

    pr = np.genfromtxt(infile, dtype='str')
    print(pr)
    print(pr[:,2].astype(float))
    max_d = np.max(pr[:,2].astype(float))
    print(np.max(pr[:,2].astype(float)))

    pr[:,2]=np.power(-np.log(pr[:,3].astype(float)), 2)
    print(pr)
    # !! incomplete

##-------------------------------------------------------------------------------------
# gene_align_*: remove abnormal AA alphabets gaps = ['.','-',' ','*']
#             : reformat fasta headers; unique ID + transcription ID
##-------------------------------------------------------------------------------------
#
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



def foo(args):
    print(args)

if __name__=='__main__':
    cp.dispatch(__name__)
