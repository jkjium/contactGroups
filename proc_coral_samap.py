import commp3 as cp
import numpy as np
import pandas as pd
from itertools import groupby
from scipy.spatial import distance
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


# h5ad assembling functions  ---------------------------------------------------------------
# append clustering assignments to an h5ad file
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

# alignment score functions --------------------------------------------------------------
# For generate heatmap
# $ ls *prost*family*|awk 'BEGIN {printf "python proc_coral_samap.py combine_scorefiles \""} {printf "%s ", $1} END{printf "\" outfile"}'
# $ ls *blast*family*|awk 'BEGIN {printf "python proc_coral_samap.py combine_scorefiles \""} {printf "%s ", $1} END{printf "\" outfile"}'
# outfile: tsv file of combined alignment scores with column ordered by cell type (family) labels
# visualize using the following R code
'''
library(pheatmap)
infile="outfile"
h=pheatmap(t(read.table(infile, row.names=1, header=TRUE, sep='\t')),cluster_rows=F)
source("save_pheatmap.R")
save_pheatmap(h, paste(infile,".png",sep=""), width=12,height=8)
'''
def combine_scorefiles(args):
    assert len(args)==2, 'Usage: python proc_coral_samap.py combine_scores_tsv "list of .tsv files" outfile'
    tsvfiles = [fn.strip() for fn in args[0].split(' ') if fn!='']
    cp._info('%s files loaded.' % len(tsvfiles))
    outfile = args[1]

    scores = pd.concat([pd.read_csv(fn, sep='\t', index_col=0).iloc[:72,72:] for fn in tsvfiles], axis=1)
    c=list(scores.columns)
    c.sort(key=lambda n: n.split('_')[1]) # sort by call type families
    scores[c].to_csv(outfile, sep='\t', index=True, header=True)
    cp._info('save all scores %s to %s' % (scores[c].shape, outfile))


# samap general pipeline functions -------------------------------------------------------
# for stage3.sandbox
# run samap for all ad vs {sp, hy, xe, nv}
def _ppl_run_samaps(map_opt='blast'):
    maps = 'maps.'+map_opt+'/'
    samaps={}

    assignments = {'ad':'Jake', 'sp':'cell_type_family'}
    sm = _run_samap_with_sam_obs('ad', '01.adig.counts_pr.cluster_info.h5ad', 'sp', '01.spis.counts_pr.cluster_info.h5ad', map_p=maps, assignments=assignments)
    save_samap(sm, '02.adsp.samap.jake_cell_type_family.'+map_opt+'.pkl')
    D,MappingTable = get_mapping_scores(sm,assignments,n_top = 0)
    MappingTable.to_csv('02.adsp.alignment_score.ad_jake.sp_cell_type_family.'+map_opt+'.tsv', sep='\t', index=True, header=True)
    samaps['adsp'] = sm


    assignments = {'ad':'Jake', 'hy':'cell_type_family'}
    sm = _run_samap_with_sam_obs('ad', '01.adig.counts_pr.cluster_info.h5ad', 'hy', '01.hvul.counts_pr.cluster_info.h5ad', map_p=maps, assignments=assignments)
    save_samap(sm, '02.adhy.samap.jake_cell_type_family.'+map_opt+'.pkl')
    D,MappingTable = get_mapping_scores(sm,assignments,n_top = 0)
    MappingTable.to_csv('02.adhy.alignment_score.ad_jake.hy_cell_type_family.'+map_opt+'.tsv', sep='\t', index=True, header=True)
    samaps['adhy'] = sm


    assignments = {'ad':'Jake', 'xe':'cell_type_family'}
    sm = _run_samap_with_sam_obs('ad', '01.adig.counts_pr.cluster_info.h5ad', 'xe', '01.xesp.counts_pr.cluster_info.h5ad', map_p=maps, assignments=assignments)
    save_samap(sm, '02.adxe.samap.jake_cell_type_family.'+map_opt+'.pkl')
    D,MappingTable = get_mapping_scores(sm,assignments,n_top = 0)
    MappingTable.to_csv('02.adxe.alignment_score.ad_jake.xe_cell_type_family.'+map_opt+'.tsv', sep='\t', index=True, header=True)
    samaps['adxe'] = sm


    assignments = {'ad':'Jake', 'nv':'cell_type_family'}
    sm = _run_samap_with_sam_obs('ad', '01.adig.counts_pr.cluster_info.h5ad', 'nv', '01.nvec.counts_pr.cluster_info.h5ad', map_p=maps, assignments=assignments)
    save_samap(sm, '02.adnv.samap.jake_cell_type_family.'+map_opt+'.pkl')
    D,MappingTable = get_mapping_scores(sm,assignments,n_top = 0)
    MappingTable.to_csv('02.adnv.alignment_score.ad_jake.nv_cell_type_family.'+map_opt+'.tsv', sep='\t', index=True, header=True)
    samaps['adnv'] = sm

    return samaps
    

# call samap main procedure
# sn1,sn2: species names
# sf1,sf2: sam object file names .h5ad (with pre-calculated clustering assignments)
# map_p: homolog graph path
# resolutions: leiden clustering resolution for each species
# assignments: = {"bf":"cluster", "mm":"tissue_celltype"}
# sm = _run_samap_with_sam_obs('ad', '00.adsp.adig.counts_pr.h5ad', 'sp', '00.adsp.spis.counts_pr.h5ad', map_p='maps')
def _run_samap_with_sam_obs(sn1, sf1, sn2, sf2, map_p='maps/', resolutions=None, assignments=None):
    from samap.mapping import SAMAP
    from samap.analysis import (get_mapping_scores, GenePairFinder, sankey_plot, chord_plot, CellTypeTriangles, ParalogSubstitutions, FunctionalEnrichment, convert_eggnog_to_homologs, GeneTriangles)
    from samap.utils import save_samap, load_samap
    from samalg import SAM
    import pandas as pd

    cp._info('loading SAM objects.')
    sam1 = SAM()
    sam1.load_data(sf1)
    sam2 = SAM()
    sam2.load_data(sf2)
    sams = {sn1: sam1, sn2: sam2}

    cp._info('running SAMap %s%s...' % (sn1,sn2))
    # keys = {"bf":"cluster", "mm":"tissue_celltype"}
    sm = SAMAP(sams, f_maps=map_p, resolutions=resolutions, keys=assignments)

    # neigh_from_keys = {'bf':True, 'mm':True}
    neigh_from_keys = dict((k,True) for k in assignments.keys()) if assignments!=None else None
    sm.run(pairwise=True, neigh_from_keys=neigh_from_keys)
    #save_samap(sm, 'ad_nv.samap.pkl')
    cp._info('SAMap %s%s done.' %(sn1,sn2))
    return sm


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



if __name__=='__main__':
    cp.dispatch(__name__)
