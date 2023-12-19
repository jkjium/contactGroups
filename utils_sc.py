import commp3 as cp
#import numpy as np
#import scipy as sp

import scanpy as sc


def leiden(args):
    sc.settings.verbosity = 0 # minimizing output
    sc.settings.set_figure_params(dpi=80, facecolor="white", frameon=False)

    adata = sc.read_h5ad("hydra.h5ad")

    sc.pp.neighbors(adata, n_pcs=30)
    sc.tl.umap(adata)

    sc.tl.leiden(adata, key_added="leiden_res0_25", resolution=0.25)
    sc.tl.leiden(adata, key_added="leiden_res0_5", resolution=0.5)
    sc.tl.leiden(adata, key_added="leiden_res1", resolution=1.0)

    sc.pl.umap(adata, color=["leiden_res0_25", "leiden_res0_5", "leiden_res1"], legend_loc="on data",)

    print(args)

def loadmatrix(args):
    assert len(args)!=4, 'Usage: python utils_sc.py loadmatrix x,y,data'
    pass

if __name__=='__main__':
    cp.dispatch(__name__)
