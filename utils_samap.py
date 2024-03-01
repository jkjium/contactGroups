import commp as cp
import numpy as np
import pandas as pd
import scipy as sp
from itertools import groupby
from scipy.spatial import distance
import matplotlib.pyplot as plt

try:
    import scanpy as sc
    from samap.mapping import SAMAP, prepare_SAMap_loadings
    from samap.analysis import (get_mapping_scores, GenePairFinder, sankey_plot, find_cluster_markers)
    from samap.utils import save_samap, load_samap
    from samalg import SAM    
    from bokeh.plotting import show
    import holoviews as hv
except ImportError:
    cp._info('ignore absent libraries')

# return euclidean distance between two vectors
def _euclidean(v1,v2):
    return np.linalg.norm(v1 - v2)

def _cosine(v1,v2):
    return np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))

# calculate the mean expression vector of a set of cells in a sam object
def _mean_expression(s, c_indices):
    # slice expression matrix by cluster member
    cluster_expression_matrix = s.adata[s.adata.obs.index.isin(c_indices)]
    return cluster_expression_matrix.X.mean(axis=0)

# calcualte wpca for sam object 
# result is stored in s.adata.varm['wpca']
# s: sam object
def _calc_varm_wpca(s):
    from sklearn.preprocessing import StandardScaler
    st = StandardScaler(with_mean=False) 

    ss = std.fit_transform(s.adata.X)
    W =s.adata.var['weights'].values
    ws = ss.multiply(W[None,:]).tocsr()
    wpca = ws.dot(s.adata.varm['PCs'])
    mu = ws.mean(0).A.flatten()

    s.adata.varm['wpca'] = wpca - mu

# must project data first project to pc space first
def _mean_expression_wpca(s, c_indices):
    pass


# s: sam object
# cluster_label: {'leiden_clusters', 'cell_type', ...}
# df = s.adata.obs['cluster_name'].copy()
# assignments = pd.DataFrame({'_neighbor_joining': df})
def _neighbor_joining(s, cluster_name, _fn_cluster_vectors=_mean_expression, _fn_metric=_euclidean): 
    assignments = pd.DataFrame({cluster_name: s.adata.obs[cluster_name].copy()})

    # get ordinal indices for constructing linkage data structure
    assignments['ordinal_index'] = pd.factorize(assignments[cluster_name])[0]

    # get mapping: m[ordinal_index]= original_label
    oi2label = dict((ordinal_i,cluster_i) for cluster_i, ordinal_i in assignments.value_counts().index)

    # get all cluster representation vectors
    cell_groups = assignments.groupby('ordinal_index').groups
    cv = dict((k, _fn_cluster_vectors(s, cell_groups[k])) for k in cell_groups)
       
    # initialize distance matrix
    n = len(cell_groups)
    distances = dict(((i,j), _fn_metric(cv[i],cv[j])) for i in range(n) for j in range(i+1,n))

    linkage=[]
    i=0
    while i<(n-1):
        # find closest pair i,j
        min_value = min(distances.values())
        c1, c2 = [k for k, v in distances.items() if v == min_value][0]
        print(i,' merging: ',c1,c2)

        # merge i,j to a new cluster
        new_label = n+i
        cell_groups[new_label] = cell_groups[c1].union(cell_groups[c2])
        linkage.append([c1,c2, distances[(c1,c2)], len(cell_groups[new_label])])

        # remove merged clusters
        d_to_remove = [k for k in distances if (c1 in k) or (c2 in k)]
        for p in d_to_remove:
            del distances[p]
        del cell_groups[c1]
        del cell_groups[c2]
        del cv[c1]
        del cv[c2]

        # update distances between existing labels and new_label
        new_cvec = _fn_cluster_vectors(s, cell_groups[new_label])
        distances.update(dict(((k, new_label), _fn_metric(new_cvec, cv[k])) for k in cv))
        cv[new_label] = new_cvec

        i+=1

    linkage_order = np.array(linkage)
    linkage_order[linkage_order[:, 0].argsort()]
    return oi2label, linkage_order

# df_table: flat homology table in pandas dataframe type
# pl.dd_Smed_v4_424_0_1 hy.t18073aep 77.273 374 85 0 135 1256 240 1361 0.0 699
# c1,c2 column pair indicating x,y ids
# c_value: homology weight {bit score}
# max_dup: if True, save only one with the maximum value when duplications found otherwise sum up all 
# reciprocate: keep only bi-directional homology 
# return scipy sparse matrix: gnnm, labels: gns, labels_by_sid: gns_dict
def _flat2spmat(df_table, c1, c2, c_value, max_dup=False, reciprocate=True):
    # pd.groupby to handle duplications
    nptable = np.array(df_table.loc[df_table.groupby([c1,c2])[c_value].idxmax()]) if max_dup==True else np.array(df_table)
    # get ordered unique labels (genes) then convert genes to integers
    labels = np.sort(np.unique(np.concatenate((nptable[:,c1], nptable[:,c2]))))
    ix = labels.searchsorted(nptable[:,c1]) 
    iy = labels.searchsorted(nptable[:,c2])
    #!!! use pd.factorize() instead

    # gnnm is asymetrical
    gnnm = sp.sparse.coo_matrix((nptable[:,c_value].astype(float), (ix,iy)), shape=(labels.size, labels.size)).tocsr()
    gnnms = (gnnm + gnnm.T) / 2
    if reciprocate:
        gnnm.data[:]=1
        gnnms = gnnms.multiply(gnnm).multiply(gnnm.T).tocsr()

    # get gns_dict: separate labels by sid
    sids = np.array([s.split('_')[0] for s in labels]) 
    sep_idx = np.where(sids[:-1]!=sids[1:])[0][0]+1 # find where sid changes
    gns_dict={sids[sep_idx-1]: labels[:sep_idx], sids[sep_idx]: labels[sep_idx:]}

    return gnnms, labels, gns_dict

# thresholding scipy sparse matrix 
# gnnm: scipy.sparse.coo_matrix; thr = scaler
# return: symmetrical filtered gnnm
def _filter_spmat(gnnm, thr=0.25):
	# threshold by row max
    thresholds = gnnm.max(axis=1).A.flatten() * thr
    x, _ = gnnm.nonzero()
    gnnm.data[gnnm.data < thresholds[x]] = 0 
    gnnm.eliminate_zeros()
	
	# patching values by letting (i,j)=(j,i)!=0 
    x,y = gnnm.nonzero()
    z = gnnm.data
    gnnm = gnnm.tolil()
    gnnm[y,x] = z
    gnnm.tocsr()
    return gnnm # symmetrical


# generate sam objects from anndata file
def _generate_sam(h5file, leiden_res=3, outfile=None):
    preprocessing= 'StandardScaler'
    weight_PCs = False

    # run SAM ppl
    sam=SAM()
    sam.load_data(h5file)
    sam.preprocess_data(
        sum_norm="cell_median",
        norm="log",
        thresh_low=-1.0,
        thresh_high=-1.96,
        min_expression=0,
    )
    sam.run(
        preprocessing=preprocessing,
        npcs=99,
        weight_PCs= weight_PCs,
        k=19,
        n_genes=2999,
        weight_mode='rms'
    )

    # calculating leiden clustering
    sam.leiden_clustering(leiden_res)

    # calculate A: weighed PCA
    A, _ = sam.calculate_nnm(
        n_genes=sam.adata.shape[0],
        preprocessing=preprocessing,
        npcs=299,
        weight_PCs=weight_PCs,
        sparse_pca=True,
        update_manifold=False,
        weight_mode='dispersion'
    )
    sam.adata.varm["PCs_SAMap"] = A

    if outfile is not None: 
        sam.save_anndata(outfile)
        cp._info('save sam h4ad to %s' % outfile)


# procedures ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def ut_flat2spmat(args):
    # unit_test/all_maps.txt
    nptable = np.loadtxt(args[0],dtype='str')

    gnnm_obj = _flat2spmat(nptable, 0,1,11,' ')
    gns = gnnm_obj['gns']
    gns_dict = dict((s[0],list(s[1])) for s in groupby(gns, lambda x: x.split('.')[0]))
    gnnms = gnnm_obj['gnnm']

    np.savetxt('mygnnm.txt', gnnms.A, fmt='%.3f')
    np.savetxt('mygnnm.tick.txt', gns, fmt='%s')
    cp._info('To generate a heatmap: python utils_vis_sm.py heatmap mygnnm.txt mygnnm.tick.txt mygnnm.tick.txt')

# append clustering assignments to an h5ad file
# nan will be renamed by unclustered_label
def _append_assignments(h5, clusterfile, cluster_colnames=['cell_type'], index_name='index', unclustered_label='unclustered'):
    cluster_assignments = pd.read_csv(clusterfile, header=None, sep='\t', index_col=0, names=cluster_colnames)
    cluster_assignments.index.name = index_name 
    h5.obs=h5.obs.merge(cluster_assignments, how='left', left_index=True, right_index=True)    
    for c in cluster_colnames:
        h5.obs[c].fillna(unclustered_label, inplace=True)

_pipeline_samap = {
    # prepare sam objects
    'prepare_input_sams': _generate_sam,
    # optional 
    'append_cluster_assignments': _append_assignments,
    # homology hits from blast or prost
    'generate_homology_graph': _flat2spmat,
    # remove weak homologous gene connections 
    'threshold_homology_graph': _filter_spmat,

}

if __name__=='__main__':
    cp.dispatch(__name__)