import numpy as np
import commp as cp


# implementation of "eigenvector centrality for characterization of protein allosteric pathways" ec analysis
# equation 9
# A_ij = 0 if i=j else rMI[X_i, X_j] * exp(-dij/lambda)
# where lambda is the locality controlling parameter
# input: 
# 1. ccmatfile, correlation matrix in np.loadtxt format
# 2. dmatfile, spatial distance matrix in np.loadtxt() format, 'na' means no padding, equivalent to set loc = 'inf'
# 3. loc, locality parameter
# 4. threshold, for all cc values < threshold, cc = 0
# 4. output file prefix, output three files {adjmatrix, adjmat.eigenvalue, adjmat.eigenvectors}
#def ccmat2adjmat(args):
def eca_ccmat2adjmat(args):
    assert len(args) == 5, 'Usage: python utils_graph.py ccmat2adjmat ccmatfile dmatfile loc thresholding outmatfile'
    ccmatfile = args[0]
    dmatfile = args[1]
    loc = float(args[2]) # 'inf' is allowed, means no damping effect
    threshold = float(args[3]) # threshold cutoff to set correlation value to zero
    outprefix = args[4]

    # for dynamic correlations
    #ccmat = np.abs(np.loadtxt(ccmatfile))
    ccmat = np.loadtxt(ccmatfile)

    #rccmat = np.loadtxt(ccmatfile)
    #ccmat = (rccmat - rccmat.min()/(rccmat.max() - rccmat.min()))
    

    # inf means no thresholding
    if threshold != float('inf'):
        cp._info('thresholding with %.4f cutoff' % threshold)
        ccmat[ccmat < threshold] = 0.0
        print ccmat
    
    # no locality damping
    if dmatfile != 'na': 
        dmat = np.loadtxt(dmatfile)
        assert ccmat.shape == dmat.shape, 'Error. dimension mismatch. ccmat: %s, dmat: %s' % (repr(ccmat.shape), repr(dmat.shape))
        dampmat = np.exp(-dmat/loc)
        adjmat = np.multiply(ccmat, dampmat)
        cp._info('apply damping factor with lambda: %.4f' % loc)
    else:
        adjmat = ccmat

    # set A{i=j} = 0.0
    np.fill_diagonal(adjmat, 0.0)
    
    # eigen decomposition    
    e, v = np.linalg.eig(adjmat)
    v = v[:, e.argsort()]
    e.sort()

    # output result adjmatrix, eigenvalues, eigenvectors, centrality (eigenvector corresponding to the largest eigenvalue)
    np.savetxt(outprefix+'.adjmat', adjmat, fmt = '%.4f')
    np.savetxt(outprefix+'.adjmat.e', e, fmt = '%.4f')
    np.savetxt(outprefix+'.adjmat.v', v, fmt = '%.4f')
    np.savetxt(outprefix+'.adjmat.ec', v[:,-1], fmt = '%.4f')
    cp._info('save output to %s {.adjmat, .adjmat.e, .adjmat.v, .adjmat.ec}' % outprefix)


# thresholding a correlation matrix into a adjacency matrix
# ccmatfile: input correlation square matrix
# thresholding: cutoff
# abs: absolute flag: 1: take abs() to the input matrix, 0: ignore
def ccmat2adjmat(args):
    assert len(args) == 4, 'Usage: python utils_graph.py ccmat2adjmat ccmatfile thresholding abs={0,1} outadjmatfile'
    ccmatfile = args[0]
    threshold = float(args[1]) # threshold cutoff to set correlation value to zero
    absflag = int(args[2])
    outfile = args[3]

    # remove negative values?
    if absflag == 1:
        ccmat = np.abs(np.loadtxt(ccmatfile))
    else:
        ccmat = np.loadtxt(ccmatfile)

    # inf means no thresholding
    if threshold != float('inf'):
        cp._info('thresholding with %.4f cutoff' % threshold)
        ccmat[ccmat < threshold] = 0.0
        ccmat[ccmat >= threshold] = 1.0
        #print ccmat

    adjmat = ccmat
    
    np.savetxt(outfile, adjmat, fmt = '%.4f')
    cp._info('save adjmat to %s' % outfile) 


def foo(args):
    print(args)
    pass

if __name__ == '__main__':
    cp.dispatch(__name__)