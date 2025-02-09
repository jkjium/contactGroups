import commp as cp
import numpy as np
import pandas as pd
import scipy as sp
import difflib
from itertools import groupby, combinations
from scipy.spatial import distance
from collections import OrderedDict, defaultdict
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

############################################################################################
# fucntions for orthofinder --------------------------------------------------
############################################################################################
# ortho one-to-one
# extract one-to-one orthologs from orthofinder outputs
# kjia@DESKTOP-K0I3VAM ~/workspace/foxd/orthofinder/Orthologues_ad 2024-12-15 18:31:00
# $ head ad__v__nt.tsv
# Orthogroup      ad      nt
# OG0000000       adig-s0038.g19  SVEP1-like-40
# OG0000000       adig-s0006.g329, adig-s0092.g75, adig-s0037.g103
# !! second column (ad) can have duplications, need to merge into a 1toMany list
# output: dict[ad] = nt
def _ortho_1to1(tuple_list):
    # merge potential 1toMany terms
    dict_1tomany = defaultdict(list)
    for t in tuple_list:
        query = t[1].split(',')
        hit = t[2].split(',')
        if(len(query)==1 and len(hit)==1):
            dict_1tomany[query[0]].append(hit[0])
    # remove 1 to many 
    return dict ((q, dict_1tomany[q][0]) for q in dict_1tomany if len(dict_1tomany[q])==1)

# load all orthofinder files and output 
def ortho_1to1_all(args):
    assert len(args) == 2, 'Usgae: python proc_foxd.py ortho_1to1_all stubfile outprefix'
    stubfile = args[0]
    outpref = args[1]
    ds = {}
    ds['names'] = cp.loadlines(stubfile)
    for n in ds['names']:
        ds[n] = _ortho_1to1(cp.loadtuples(n, '\t')[1:])
    # find common keys
    common_set = set(ds[ds['names'][0]].keys())
    for n in ds['names'][1:]:
        common_set &= set(ds[n].keys())
    cp._info('%d common keys found' % len(common_set))
    for n in ds['names']:
        outfile = '%s.%s.tsv' % (outpref, n) 
        with open(outfile , 'w') as fout:
            fout.write('%s\n' % '\n'.join(['%s\t%s' % (k, ds[n][k]) for k in common_set]))
        cp._info('save filtered 1to1 ortholog to %s' % outfile)

# ut _ortho_1to1()
def ortho_1to1(args):
    assert len(args) == 2, 'Usgae: python proc_foxd.py ortho_1to1 ad__v__nt.tsv outfile'
    infile = args[0]
    outfile = args[1]
    ortho_pair_list = cp.loadtuples(infile, '\t')[1:]
    out_1to1_dict = _ortho_1to1(ortho_pair_list)
    with open(outfile, 'w') as fout:
        fout.write('%s\n' % '\n'.join(['%s\t%s' % (k, out_1to1_dict[k]) for k in out_1to1_dict]))
    cp._info('save %d 1to1 ortho pair to %s' % (len(out_1to1_dict), outfile))

# filter the shared 1to1 orthologs from orthofinder
# some genes may not exist in seurat objects
def filter_1to1(args):
    assert len(args) == 3, 'Usage: python proc_coral_samap.py filter_1to1 gene_names.stub ortho_1to1.tsv ortho_1to1.filtered.tsv'
    stubfile = args[0] # gene names from seruat objects
    orthofile = args[1] # shared 1to1 orthologs from orthofinder
    outfile = args[2] 

    # load all seurat object genes into a list
    genes = []
    for f in cp.loadlines(stubfile):
        genes.extend(cp.loadlines(f))
    cp._info('load %d genes' % len(genes))

    # check 1to1 ortholog genes line by line
    outlist = []
    for ortho_row in cp.loadtuples(orthofile, '\t'):
        if len([g for g in ortho_row if g in genes]) == len(ortho_row):
            outlist.append('\t'.join(ortho_row))
    
    with open(outfile, 'w') as fout:
        fout.write('%s\n' % '\n'.join(outlist))
    cp._info('save %d filtered 1to1 orthologs to %s' % (len(outlist), outfile))


############################################################################################
# fucntions for visualization --------------------------------------------------
############################################################################################
# df: three column df, c1: name, c2: value, c3: value
# for plotting sym/apo ratio
# df = pd.read_csv('info.gastrodermis.txt', sep=' ', header=None, names=['clusterID', 'apo', 'sym'])
def pie_charts(df, fs=14):
    import matplotlib.pyplot as plt
    n = len(df)
    fig, axs = plt.subplots(1, n, figsize=(n * 4, 4))
    for index, row in df.iterrows():
        sizes = [row['apo'], row['sym']]
        labels = ['apo', 'sym']
        axs[index].pie(sizes, labels=labels, autopct='%1.1f%%', startangle=90, textprops={'fontsize': fs})
        axs[index].set_title(row['clusterID'], fontsize=fs)
        axs[index].axis('equal')
        #axs[index].text(0, -1.2, row['clusterID'], fontsize=12, ha='center')
    plt.tight_layout()
    plt.show()


def extract_seq_by_list(args):
    assert len(args) == 3, 'Usage: python proc_coral_samap.py extract_seq_from_list seq.fasta names.txt outfile'
    from tqdm import tqdm
    infile = args[0]
    stubfile = args[1]
    outfile = args[2]

    stub = cp.loadlines(stubfile)
    c=0
    outlist=[]
    for s in cp.fasta_iter(infile):
        sarr = s[0].split('.')
        gname= '.'.join(sarr[:2])
        #print(gname)
        if gname in stub:
            outlist.append('>%s\n%s' % (s[0], s[1]))
            c+=1
            print(c)
    with open(outfile, 'w') as fout:
        fout.write('%s\n' % '\n'.join(outlist))
    cp._info('save sequences to %s' % outfile)


############################################################################################
# fucntions for gene annotations --------------------------------------------------
############################################################################################
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
    cp._info('save to %s' % outfile)

# for only gene name extraction
def get_gene_name(args):
    assert len(args)==2, 'Usage: python proc_coral_samap.py get_gene_name curl.json.txt outfile'
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
            gid = jo['results'][0].get('genes', '_NA_')
            gn = gid[0].get('geneName', '_NA_') if gid!='_NA_' else '_NA_'
            gnv = gn.get('value', '_NA_') if gn!='_NA_' else '_NA_'
        except Exception as e:
            print('exception: [%s]' % e)
            print(jsonline)
        outstr= '%s\t%s' % (uid, gnv)
        print('extract: %s' % outstr)
        outlist.append(outstr)
        jsonline=''
    fp.close()
    with open(outfile, 'w') as fout:
        fout.write('%s\n' % '\n'.join(outlist))
    cp._info('save to %s' % outfile)




# extract human, drome, top1 from *allhits.cell.uid.pn.gd.gn.tax.tn.dis.tsv
# gene names are sorted by prost distance
def ann_merge(args):
    assert len(args)==2, 'Usage: python proc_coral_samap.py ann_merge in.allhits.cell.uid.pn.gd.gn.tax.tn.dis.tsv outfile'
    infile = args[0]
    outfile = args[1]
    # load data
    df = pd.read_csv(infile, sep='\t', names=['Query_gene_ID','Uniprot_ID','Protein_name','Gene_name','Taxon_ID','Taxon_name','PROST_distance'])
    # filter _NA_ gene names
    df=df.loc[df['Gene_name']!='_NA_']

    #f=lambda x: ','.join(list(x))
    f = lambda x: ','.join(x.sort_values(by='PROST_distance', key=lambda col: col.astype(float))['Gene_name'].tolist())

    # extract human hits
    dh = df.loc[df['Taxon_ID']== 9606]
    #dh = dh.groupby('Query_gene_ID')['Gene_name'].apply(f).reset_index()
    dh = dh.groupby('Query_gene_ID').apply(f).reset_index()
    dh.set_index('Query_gene_ID', inplace=True)
    #dh.rename(columns={'Gene_name': 'PROST_gene_name(Human)'}, inplace=True)
    dh.rename(columns={0: 'PROST_gene_name(Human)'}, inplace=True)

    # extract drome hits
    dd = df.loc[df['Taxon_ID']== 7227]
    #dd = dd.groupby('Query_gene_ID')['Gene_name'].apply(f).reset_index()
    dd = dd.groupby('Query_gene_ID').apply(f).reset_index()
    dd.set_index('Query_gene_ID', inplace=True)
    #dd.rename(columns={'Gene_name': 'PROST_gene_name(Drome)'}, inplace=True)
    dd.rename(columns={0: 'PROST_gene_name(Drome)'}, inplace=True)

    # extract gene name, tax name from the top1 hit
    #da = df.groupby('Query_gene_ID').apply(lambda x: sorted(list(zip(x['Uniprot_ID'], x['Protein_name'], x['Gene_name'], x['Taxon_ID'], x['Taxon_name'], x['PROST_distance'])), key=lambda tup: tup[5])[0])
    da = df.groupby('Query_gene_ID').apply(lambda x: sorted(list(zip(x['Query_gene_ID'], x['Uniprot_ID'], x['Protein_name'], x['Gene_name'], x['Taxon_ID'], x['Taxon_name'], x['PROST_distance'])), key=lambda tup: tup[6])[0])
    da = da.reset_index(name='PROST_gene_name(Top1)')

    # output top1 flat tsv
    dt = pd.DataFrame(da['PROST_gene_name(Top1)'].tolist(), columns=['Query_gene_ID', 'Uniprot_ID','Protein_name','Gene_name','Taxon_ID','Taxon_name','PROST_distance'])
    dt.set_index('Query_gene_ID', inplace=True)
    dt.to_csv('out.top1.tsv', sep='\t', index=True, header=True)
    cp._info('save top1 df to out.top1.tsv')

    da['PROST_gene_name(Top1)'] = da['PROST_gene_name(Top1)'].apply(lambda x: '%s, %s' % (x[3], x[5]))
    da.set_index('Query_gene_ID', inplace=True)
    # merge dh,dd,da
    outdf = da.merge(dh, how='left', left_index=True, right_index=True).merge(dd, how='left', left_index=True, right_index=True) 
    outdf.fillna('_NA_', inplace=True)
    outdf.to_csv(outfile, sep='\t', index=True, header=True)
    cp._info('save merged df to %s' % outfile)

    # select names for each query gene ID
    # Truncate to 50 characters
    def _choose_alias(row):
        v1 = row['PROST_gene_name(Human)']
        v2 = row['PROST_gene_name(Drome)']
        v3 = row['PROST_gene_name(Top1)']
        
        _t = lambda x: x if len(x)<=50 else '%s..' % x[:50]

        if v1 != "_NA_":
            return 'hs-%s' % _t(v1)
        elif v2 != "_NA_":
            return 'dm-%s' % _t(v2)
        else:
            return _t(v3)

    sele_df = pd.DataFrame({'prost_alias': outdf.apply(_choose_alias, axis=1)})
    sele_df.index = outdf.index
    sele_df.index.names=['geneID'] # make consistent with emapp final results for merging
    outfile = '21.geneID.geneAlias.tsv'
    sele_df.to_csv(outfile, sep='\t', index=True, header=True)
    cp._info('save truncated(50) alias df to %s' % outfile)


# merge prost alias into emapper final names
def ann_merge_alias(args):
    assert len(args) == 3, 'Usage: python proc_coral_samap.py ann_merge_alias aten_emapper_names_table_final.tsv 21.geneID.geneAlias.tsv aten_gene_alias.tsv'
    emapper_file = args[0]
    prost_file = args[1]
    outfile = args[2]

    emapper_df = pd.read_csv(emapper_file, sep='\t')
    emapper_dict = emapper_df.set_index('geneID')['final_emapper_name'].to_dict()

    prost_df = pd.read_csv(prost_file, sep='\t')
    prost_dict = prost_df.set_index('geneID')['prost_alias'].to_dict()

    out_dict = {}
    # make sure all keys are "_" seperated as "adig_s0001.g1"
    key_stub = set(emapper_dict.keys()).union(set(prost_dict.keys()))
    for k in key_stub:
        if k in emapper_dict and k in prost_dict:
            if emapper_dict[k] == k.replace('_', '-'): # if emapper has no hit, append prost result
                out_dict[k] = '%s %s' % (emapper_dict[k], prost_dict[k])
            else:
                out_dict[k] = emapper_dict[k]
        elif k in emapper_dict: # only has emapper hit
            out_dict[k] = emapper_dict[k]
        elif k in prost_dict: # only has prost hit
            out_dict[k] = prost_dict[k]

    with open(outfile, 'w') as fout:
        fout.write('geneID\tgene_alias\n')
        fout.write('%s\n' % '\n'.join(['%s\t%s' % (k, out_dict[k]) for k in out_dict]))
    cp._info('save merged gene names to %s' % outfile)



############################################################################################
# fucntions for homology graph comparison --------------------------------------------------
############################################################################################
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
def _samh5ad(h5file,outfile, umap=True, leiden=True):
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
        weight_mode='rms',
        umap_flag = umap
    )
    if leiden == True:
        print('Calculating leiden clusters ...')
        sam.leiden_clustering(res=3)
    else:
        print('skip leiden clustering')
    
    if "PCs_SAMap" not in sam.adata.varm.keys():
        print('Calculating PCA ...')
        prepare_SAMap_loadings(sam)  

    sam.save_anndata(outfile)
    cp._info('save sam h5ad to %s' % outfile)

# slice adata obj by labels of cluster_name in obs 
# all columns in var will be removed
# h5: input anndata object
# target_obs_column:'Jake'
# labels = ['Jake_10', 'Jake_17', 'Jake_20', 'Jake_22', 'Jake_27', 'Jake_46', 'Jake_54']; corresponds to xe_gastrodermis cell type family
# return: new anndata object contains expression counts for selected cells
def _slice_by_cluster(h5, target_obs_column, labels, obs_droplist=None):
    c = h5.obs[target_obs_column]
    filters = c[c.isin(labels)].index 
    sub_h5 = h5[filters]
    obs = sub_h5.obs.drop(columns=obs_droplist) if obs_droplist is not None else sub_h5.obs
    var = sub_h5.var.drop(columns=sub_h5.var.columns) if sub_h5.var.columns!=0 else sub_h5.var
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
    #cluster_assignments = pd.read_csv(clusterfile, header=None, sep=delimiter, index_col=0, names=cluster_colnames)
    cluster_assignments = pd.read_csv(clusterfile, sep=delimiter, index_col=0)
    cluster_assignments.columns = cluster_colnames
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
    cp._info('loading mtx ...')
    h5 = _mtx2h5ad(args[0], args[1], args[2])
    cp._info('appending clustering assignments ...')
    _append_assignments(h5, args[3], [args[4]])
    cp._info('done.')
    h5.write_h5ad(args[5])
    cp._info('save to %s' % args[5])



############################################################################################
# alignment score functions --------------------------------------------------------------
############################################################################################

# find the longest substring between two strings
def _sharedstr(s1, s2):
    m=difflib.SequenceMatcher(None, s1, s2).find_longest_match(0, len(s1), 0, len(s2))
    return s1[m.a:m.a+m.size]

# find shortest common substring by comparing the first one with others for the current set 
def _min_comm_str(sset, cutoff=4):
    ms = min([_sharedstr(sset[0], si) for si in sset[1:]], key=len)
    return ms if len(ms)>=cutoff else ''

# ss: a list of strings
# num_summary: integer cutoff for top summaries
def _get_summaries(ss, num_summary=2):
    # iterate all rth order subsets
    # find summary for the each order
    for r_order_set, ord in [(combinations(ss, r), r) for r in range(len(ss), 1, -1)]: 
        # get {num_summary} longest common patterns among the current rth order subsets
        comm_strs = sorted([_min_comm_str(tuple) for tuple in r_order_set], key=len, reverse=True)[:num_summary]
        # remove empty strings
        ret_strs = [s for s in comm_strs if len(s)!=0]
        if len(ret_strs)==0:
            continue
        elif len(ret_strs)==1:
            #return '(%d) - %s' % (ord, comm_strs[0])
            return '%s (%d)' % (comm_strs[0], ord)
        else:
            #return '(%d) - %s' % (ord, ' - '.join(comm_strs))
            return '%s (%d)' % (' - '.join(comm_strs), ord)
    return ''


# for generating the cluster name alias in dotplots
# score file: combined scores, output from
# cutoff: alignment score cutoff
# trim_str: scorefile contains name like: "ad_g1.sym", need to trim "ad_" to match wgcna dotplot x-axis labels
# python proc_coral_samap.py cluster_alias_with_cutoff 04.heatmap.at2all.merged_clusters_cell_type_family.blast.tsv 0.4 at_ 05.at2all_merged_clusters{_alias_04.tsv/_summary_04.tsv}
def cluster_alias_with_cutoff(args):
    assert len(args) == 4, 'Usage: python proc_coral_samap.py cluster_alias_with_cutoff combined_score.tsv cutoff ad_ outprefix'
    scorefile=args[0]
    cutoff=float(args[1])
    trim_str = args[2]
    outprefix = args[3]

    m = pd.read_csv(scorefile, sep='\t', index_col=0)
    result_dict = OrderedDict()
    summary_dict = OrderedDict()
    for index, row in m.iterrows():
        row_sum = row.sum()
        if row_sum == 0:
            result_dict[index] = [] 
        else:
            #normalized_row = row / row_sum
            normalized_row = row # for future normalization
            # col_score = [(col.name, score),...]
            col_score = list(normalized_row.items())
            # Sort column-score pairs by scores 
            sorted_col_score = sorted(col_score, key=lambda x: x[1], reverse=False)
            mlist = [(n, v) for n, v in sorted_col_score if v >= cutoff]
            result_dict[index] = mlist
            # get cluster assignment names
            # 'ad_neuro_gland' -> 'neural_gland'
            ss = ['_'.join(t[0].split('_')[1:]) for t in mlist]
            summary_dict[index] =  _get_summaries(ss)
            #print(index, result_dict)
    
    # output raw samap alias
    outfile = '%s_alias_%.2f.tsv' % (outprefix, cutoff)
    with open(outfile, 'w') as fout:
        fout.write('idx\talias\n')
        for k in result_dict:
            outname = k.replace(trim_str, "")
            if len(result_dict[k]) == 0:
                fout.write('%s\t%s\n' % (outname, outname))
            else:
                alias = ''.join(['%s(%.2f) - ' % (t[0], t[1]) for t in result_dict[k]])
                fout.write('%s\t%s%s\n' % (outname, alias, outname))
    cp._info('save alias to %s' % outfile)

    # output summary alias
    outfile = '%s_summary_%.2f.tsv' % (outprefix, cutoff)
    with open(outfile, 'w') as fout:
        fout.write('idx\talias\n')
        for k in result_dict:
            outname = k.replace(trim_str, "")
            if len(result_dict[k]) == 0:
                fout.write('%s\t%s\n' % (outname, outname))
            elif len(summary_dict[k])!=0:
                fout.write('%s\t%s - %s\n' % (outname, summary_dict[k], outname))
            else:
                alias = ''.join(['%s(%.2f) - ' % (t[0], t[1]) for t in result_dict[k]])
                fout.write('%s\t%s%s\n' % (outname, alias, outname))
    cp._info('save summary to %s' % outfile)



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

# python proc_coral_samap.py reorder_mat 04.sl2all.alignment_scores.blast.cell_type_sele_cell_type_family.tsv heatmap.x.sci.stub na out
# reorder the combined alignment scores to generate heatmap
# pheatmap(t(read.table(infile, row.names=1, header=TRUE, sep='\t')), cluster_rows=F, cluster_cols=F)
def reorder_mat(args):
    assert len(args)==4, 'Usage: python proc_coral_samap.py reorder_mat mat_with_label.tsv x_order.stub y_order.stub outfile'

    # row: x-axis; column: y-axis in the heatmap
    scores = pd.read_csv(args[0], sep='\t', index_col=0, header=0)
    x_labels = cp.loadlines(args[1]) if args[1]!= 'na' else list(scores.index)
    y_labels = cp.loadlines(args[2]) if args[2]!= 'na' else list(scores.columns)
    scores_ro = scores.reindex(index=x_labels, columns=y_labels)

    scores_ro.to_csv(args[3], sep='\t', index=True, header=True)
    cp._info('save re-ordered scores to %s' % args[3])



############################################################################################
# samap general pipeline functions -------------------------------------------------------
############################################################################################



# si: speicies identifier {'at','hy','xe' ...}
# h5: SAM h5ad file name
# ci: cluster assignment identifier {'cell_type_family', 'merged_clusters' ...}
# align_thr: cutoff for sankey plot
# output: score_matrix.tsv, samap_obj.pkl
# example: sm = std_samap_with_assignments('ad', '01.ad3456.sam.h5ad', 'merged_clusters', 'hy', '01.hy.sam.cluster_info3.h5ad', 'cell_type_family')
# output example: 02.adhy.samap.blast.merged_clusters_cell_type_family.pkl, 03.adhy.alignment_scores.blast.merged_clusters_cell_type_family.tsv
def std_samap_with_assignments(si1, h51, ci1, si2, h52, ci2, align_thr=0.0):
    assignments = {si1:ci1, si2:ci2}
    sm = _run_samap_with_h5(si1, h51, si2, h52, samobj=True, map_p='maps.blast/', gnnm=None, assignments=assignments, res_dk=None, gnnm_transform=True)
    save_samap(sm, '02.%s%s.samap.blast.%s_%s.pkl' % (si1,si2,ci1,ci2))
    MappingTable = samap_alignment_and_sankey(sm, si1, ci1, si2, ci2, sankey_cutoff=align_thr)
    return sm, MappingTable


# sm: samap object
# si1, si2: species IDs
# ci1, ci2: cluster assignment IDs
# sankey_cutoff
# output: alignment_scores.tsv used for heatmap
# return: alignment score matrix
# pheatmap
# len(np.unique(sm.sams['ad'].adata.obs['merged_clusters'])) = 65
# $ l --color=none 03.ad*alignment*.blast.*cell_type_family*.tsv |awk 'BEGIN {printf "python proc_coral_samap.py combine_scorefiles \""} {printf "%s ", $1} END{printf "\" 65 04.heatmap.ad2all.merged_clusters_cell_type_family.blast.tsv"}' |sh
# 
# library(pheatmap)
# infile="04.heatmap.ad2all.merged_clusters_cell_type_family.blast.tsv"
# h=pheatmap(t(read.table(infile, row.names=1, header=TRUE, sep='\t')),cluster_rows=T)
# the order of si1,ci1,si2,ci2 must be consistent with the order when running samap
def samap_alignment_and_sankey(sm, si1, ci1, si2, ci2, sankey_cutoff=0.0):
    D,MappingTable = get_mapping_scores(sm, {si1:ci1, si2:ci2} ,n_top = 0)
    MappingTable.to_csv('03.%s%s.alignment_scores.blast.%s_%s.tsv' % (si1,si2,ci1,ci2), sep='\t', index=True, header=True)
    cp._info('save full score matrix to: 03.%s%s.alignment_scores.blast.%s_%s.tsv' % (si1,si2,ci1,ci2))
    # remove redundant info from M
    # full MappingTable will have problem when all values in a row are identical
    p=len(np.unique(sm.sams[si1].adata.obs[ci1]))
    mr = MappingTable.iloc[:p,p:]

    cmap = _sankey_cmap(si2, sm.sams[si2].adata.obs[ci2], si1, sm.sams[si1].adata.obs[ci1], M=mr, cmap_id=None, color_scheme=_cscheme_seurat_dimplot)
    show(hv.render(sankey_plot(MappingTable, align_thr=sankey_cutoff, species_order = [si1,si2], cmap=cmap)))
    return MappingTable

# call samap main procedure
# sn1,sn2: species names
# sf1,sf2: sam object file names .h5ad (with pre-calculated clustering assignments)
# map_p: homolog graph path
# resolutions: leiden clustering resolution for each species
# assignments: = {"bf":"cluster", "mm":"tissue_celltype"}
# _run_samap_with_h5('ad', '00.adig.counts.Jake_xe_gastodermis.h5ad', 'xe', '00.xesp.counts.Jake_xe_gastodermis.h5ad', samobj=False, map_p='maps.blast/', assignments=assignments)
def _run_samap_with_h5(sn1, fn1, sn2, fn2, samobj=False, map_p='maps/', gnnm=None, resolutions=None, assignments=None, res_dk=1.0, gnnm_transform=True, umap=False, niter=3):
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
    sm.run(pairwise=True, neigh_from_keys=neigh_from_keys, umap=umap, NUMITERS=niter) # umap is false for debugging
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

# input: blast format maps
#  1.  qseqid      query or source (gene) sequence id
#  2.  sseqid      subject or target (reference genome) sequence id
#  3.  pident      percentage of identical positions
#  4.  length      alignment length (sequence overlap)
#  5.  mismatch    number of mismatches
#  6.  gapopen     number of gap openings
#  7.  qstart      start of alignment in query
#  8.  qend        end of alignment in query
#  9.  sstart      start of alignment in subject
# 10.  send        end of alignment in subject
# 11.  evalue      expect value  # get minimum
# 12.  bitscore    bit score  # need normalize

# use the first two column as key to merge scores with same query,match sequence
# output: non-redundant sskey,bitscores with re-scaled bitscore
def _congruent_blast_map(mfile, remove_redundancy=True):
    df = pd.read_csv(mfile,sep='\t', names=['qseqid', 'sseqid', 'pident', 'length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore'])
    df['sskey'] = df['qseqid']+'\t'+df['sseqid']
    if remove_redundancy:
        idx = df.groupby(['sskey'])['evalue'].idxmin()
        df=df.loc[idx]
    df['bitscore']=df['bitscore']/df['bitscore'].max() # scale by maximum
    # remove unused columns
    df.drop(columns=['qseqid', 'sseqid', 'pident', 'length','mismatch','gapopen','qstart','qend','sstart','send'], inplace=True)
    df.set_index('sskey', inplace=True)
    return df


# merge blast, prost hits into one
def congruent_map(args):
    assert len(args) == 3, 'Usage: python proc_coral_samap.py congruent_map mapfile nr{0|1} outfile'
    mfile = args[0]
    nr = int(args[1])
    outfile = args[2]

    df = _congruent_blast_map(mfile, nr)
    import csv
    df.to_csv(outfile, sep='\t', index=True, header=False, quoting=csv.QUOTE_NONE, doublequote=False)
    #df.to_csv(outfile, sep='\t', index=True, header=False, doublequote=False)
    cp._info('save conguent map to %s' % outfile)


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


# aten remove the last .t1 and find the longest
def _at_parser(slist, arg='aten.pep.t1.fasta'):
    _fn_h = lambda sarr: '%s.%s' % (sarr[0], sarr[1])
    seqs = [[_fn_h(s[0].split('.')), s[0], s[1], len(s[1])] for s in slist]
    df = pd.DataFrame(seqs, columns=['output_id','old_id', 'seq', 'length'])
    df.to_csv(arg+'.stat.list', columns=['output_id', 'old_id', 'length'], sep='\t', header=False, index=False)
    idx = df.groupby(['output_id'])['length'].idxmax() # get the longest
    return df.loc[idx]


# for adig
# s[0]: adig-s0001.g1.t2 # remove the last .t2
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


# fa header:
# >TriadP8101 pep scaffold:ASM15027v1:scaffold_41:90438:90993:1 gene:TriadG8101 transcript:TriadT8101 gene_biotype:protein_coding transcript_biotype:protein_coding gene_symbol:SYT15 description:Synaptotagmin
# 15 [Source:UniProtKB/TrEMBL;Acc:B3SDG8]
#
# c1 is the same as c10 
# TriadP55444 TriadG55444
#
# kjia@DESKTOP-L0MMU09 ~/workspace/samap_general/data/tr/stage.clean 2024-03-12 10:54:21
# $ awk -F '_' '{print $1}' tr_genes.tsv |sort|uniq -c
#      1 DLX
#      1 RPS5
#      1 SECP1
#      1 SYT15
#      1 SYTNTM
#      1 SYTWLL
#  16380 Tadh
# >>> h5.X.shape: (13236, 16386)
# 16380 (Tadh_TriadG51255) + manual convert 6 (convert gene_symbol to gene)
def _tr_parser(slist, arg='Trichoplax_adhaerens.ASM15027v1.pep.all.fa'):
    symlist = ['DLX', 'RPS5', 'SECP1', 'SYT15', 'SYTNTM', 'SYTWLL']

    def _fn_get_var(h):
        sym=''
        if 'gene_symbol' in h:
            hsub = h.split('gene_symbol:')
            sym = hsub[1].split(' ')[0]
        return sym if sym in symlist else 'Tadh_TriadG%s' % h.split(' ')[0][6:]

    seqs = [[_fn_get_var(s[0]), s[0], s[1], len(s[1])] for s in slist]

    df = pd.DataFrame(seqs, columns=['output_id','old_id', 'seq', 'length'])
    df.to_csv(arg+'.stat.list', columns=['output_id', 'old_id', 'length'], sep='\t', header=False, index=False)
    idx = df.groupby(['output_id'])['length'].idxmax()
    return df.loc[idx]    


# no need to change
def _pd_parser(slist, arg='pdumv021.xloc.pep'):
    seqs = [[s[0], s[0], s[1], len(s[1])] for s in slist]
    df = pd.DataFrame(seqs, columns=['output_id','old_id', 'seq', 'length'])
    df.to_csv(arg+'.stat.list', columns=['output_id', 'old_id', 'length'], sep='\t', header=False, index=False)
    idx = df.groupby(['output_id'])['length'].idxmax()
    return df.loc[idx]

# parsing transdecoder output of pacbio transcriptome
def _pt_parser(slist, arg):
    cp._info('parsing transdecoder output of pacbio transcriptome')
    # 'BioSample_1_transcript/4901620.p1 type:internal gc:universal BioSample_1_transcript/4901620:1-228(+)'
    # s[0].split('.')[0][length("BioSample_1_transcript/"):] = 4901620
    # get only trascription ID {4901620} to groupby->get the longest
    seqs = [[s[0].split('.')[0][len("BioSample_1_transcript/"):], s[0], s[1], len(s[1])] for s in slist]
    df = pd.DataFrame(seqs, columns=['output_id','old_id', 'seq', 'length'])
    df.to_csv(arg+'.stat.list', columns=['output_id', 'old_id', 'length'], sep='\t', header=False, index=False)
    idx = df.groupby(['output_id'])['length'].idxmax()
    return df.loc[idx]

# keep the original header and get longest sequence for each transcription ID
def _default_parser(slist, arg):
    cp._info('use default parser.')
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
    assert len(args) == 4, 'Usage: python proc_coral_samap.py process_proteome proteome.fas parser_id gene_name.tsv/na outfile'

    _parser_list = { 
        'def':_default_parser, 'nt': _nt_parser, 'at': _at_parser,'ad': _ad_parser, 'sl': _sl_parser, 'xe': _xe_parser, 
        'hy': _hy_parser, 'sp': _sp_parser, 'nv': _nv_parser, 'tr': _tr_parser, 'pd': _pd_parser,
        'pacbio-transcriptome': _pt_parser

        }

    proteome_file = args[0]
    _fn_s_parser=_parser_list[args[1]]
    cp._info('Parser %s selected' % _fn_s_parser.__name__)
    outfile = args[3]

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
    df_filtered = _fn_s_parser(slist, proteome_file)
    print(df_filtered.shape, df_filtered.iloc[:5])
    #df_filtered.to_csv(proteome_file+'.filtered.list', columns=['output_id', 'old_id', 'length'], sep='\t', header=False, index=False)
    cp._info('save filtered header mapping to %s' % proteome_file+'.stat.list')

    if args[2]!='na':
        cp._info('gene sym file: %s provided. checking matches' % args[2])
        gene_syms = set(cp.loadlines(args[2]))
        fa_syms = set(df_filtered['output_id'])
        matches = gene_syms.intersection(fa_syms)

        print('\n------- gene name mapping statistics ------')
        cp._info('%d matches' % len(matches))
        missing = pd.Series([g for g in gene_syms if g not in fa_syms])
        cp._info('%d missed' % len(missing))
        if len(missing)!=0:
            missing.to_csv(args[0]+'.missing.tsv', header=False, index=False)
            cp._info('save missed gene sym to %s.missing.tsv' % args[0])
    else:
        cp._info('skip checking gene sym matches')

    # output the converted proteome
    print('------------------------------------------\n')
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
def _filter_small_clusters(df, n=2):
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
        # keys_to_remove = M.index[(M == 0).all(axis=1)].tolist() # find index with all zeros, they will cause the following problem
        # clusters2 = clusters2[~np.isin(clusters2, keys_to_remove)]

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
