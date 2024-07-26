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
def _mean_expression(s, cell_names):
    # slice expression matrix by cluster member
    cluster_expression_matrix = s.adata[s.adata.obs.index.isin(cell_names)]
    return cluster_expression_matrix.X.mean(axis=0).A.flatten()

# calcualte wpca for a non-samap sam object 
# result is stored in s.adata.varm['wpca']
# s: sam object
def _calc_varm_wpca(s):
    from sklearn.preprocessing import StandardScaler
    std = StandardScaler(with_mean=False) 

    # project normalized, weighed data on PCs
    W =s.adata.var['weights'].values
    ss = std.fit_transform(s.adata.X).multiply(W[None,:]).tocsr()
    wpca = ss.dot(s.adata.varm['PCs'])

    # project mean on PCs
    mu = ss.mean(0).A.flatten()
    mu_pca = mu.dot(s.adata.varm['PCs'])

    s.adata.uns['wpca'] = wpca - mu_pca

# must project data first project to pc space first using _calc_varm_wpca()
def _mean_expression_wpca(s, cell_names):
    n_indices = s.adata.obs.index.get_indexer(cell_names) # get position index
    if 'wpca' not in s.adata.uns:
        print('projecting data to wPC space')
        _calc_varm_wpca(s)
    sliced_wpca = s.adata.uns['wpca'][n_indices, :] # slice wpca matrix
    return sliced_wpca.mean(axis=0)

# calculate representative vector for each cluster id
# output: [[label, vector], [label2, vector], ...]
def _get_cluster_vectors(s, cluster_name, fn_cluster_vector=_mean_expression):
    return [[label, fn_cluster_vector(s, cell_names.tolist())] for label, cell_names in s.adata.obs.groupby(cluster_name, observed=False).groups.items()]

# for tree generation
# input: sam object, cluster_assignment_name, cluster_representative, pairwise_distance_measure
# output: dictionary {labels, pairwise_distance_matrix}
def _calc_distance_matrix(s, cluster_name, fn_cluster_vector=_mean_expression, metric='cosine'):
    from scipy.spatial.distance import pdist, squareform
    nptable_label_vector = np.array(_get_cluster_vectors(s, cluster_name, fn_cluster_vector), dtype=object)
    labels = nptable_label_vector[:,0]
    vectors = np.vstack(nptable_label_vector[:,1])
    distance_matrix = squareform(pdist(vectors, metric))
    return pd.DataFrame(distance_matrix, index=labels, columns=labels)


# discarded, using R::ape::nj instead
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


# np version of entropy
def _h_np(mat):
    nm = mat/(mat.sum(axis=1).reshape(-1, 1))
    lnm = np.where(nm!= 0, np.log2(nm), 0)
    return -(nm*lnm).sum(axis=1)

# scipy version of entropy
def _h_sp(mat):
    print("sp")
    nm = mat.multiply(np.reciprocal(mat.sum(axis=1))) # row sum normalized matrix
    lnm = nm.copy()
    lnm.data = np.log2(lnm.data)
    lnm.eliminate_zeros()
    return -np.asarray(nm.multiply(lnm).sum(axis=1))

# calculate expression entropy for given expression matrix
# assuming rows are cells and columns are genes
# input: .mtx matrix file
# output: single column entropy value corresponds to the order of barcodes
def expression_entropy(args):
    assert len(args) == 2, 'Usage: python utils_samap.py expression_entropy in.mtx out.vec'
    infile = args[0]
    outfile = args[1]

    mat = sp.io.mmread(infile).astype(float) 
    mat = mat.T # 10x matrix needs to be transposed
    #h_vec = _h_np(mat.A) 
    h_vec = _h_sp(mat) # if the matrix is huge
    print("examine the entropy dimension with number of cells: ")
    print(h_vec.shape)

    nfeature = np.array(mat.getnnz(axis=1)).reshape(-1,1)
    # output format: entropy row_sum
    np.savetxt(outfile, np.hstack((h_vec.reshape(-1,1), np.asarray(mat.sum(axis=1)), nfeature)), fmt='%.8f')
    cp._info('save entropy, nCounts, nFeatures to %s' % outfile)


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