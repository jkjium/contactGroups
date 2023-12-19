import commp3 as cp
import numpy as np
import scipy as sp
from itertools import groupby

# nptable: flat homology table
# pl.dd_Smed_v4_424_0_1 hy.t18073aep 77.273 374 85 0 135 1256 240 1361 0.0 699
# hy.t1aep sc.Smp_241470 34.762 210 128 5 144 755 323 529 3.50e-30 100
# c1,c2 column pair indicating x,y ids
# c_value: homology weight
# reciprocate: keep only bi-directional homology 
# return scipy sparse matrix: gnnm, labels: gns
def _flat2spmat(nptable, c1, c2, c_value, delimiter=',', reciprocate=True):
    # get ordered unique labels (genes)
    labels = np.sort(np.unique(np.concatenate((nptable[:,c1], nptable[:,c2]))))
    # convert genes to integers
    ix = labels.searchsorted(nptable[:,c1])
    iy = labels.searchsorted(nptable[:,c2])

    # gnnm is asymetrical
    gnnm = sp.sparse.coo_matrix((nptable[:,c_value].astype(float), (ix,iy)), shape=(labels.size, labels.size)).tocsr()
    gnnms = (gnnm + gnnm.T) / 2
    if reciprocate:
        gnnm.data[:]=1
        gnnms = gnnms.multiply(gnnm).multiply(gnnm.T).tocsr()
    
    gnnm_obj = {'gnnm': gnnms, 'gns':labels}
    return gnnm_obj

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


def foo(args):
    pass

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


if __name__=='__main__':
    cp.dispatch(__name__)