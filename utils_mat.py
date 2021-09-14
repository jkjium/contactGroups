import numpy as np
import commp as cp

def _trimbytick(mat, fulltick, trimtick):
    trimidlist = [fulltick.index(i) for i in (set(fulltick) - set(trimtick))]
    outmat = np.delete(mat, trimidlist, axis=0) # trim rows
    outmat = np.delete(outmat, trimidlist, axis=1) # trim columns
    return outmat


# resi list from cemat is less than dyncc mat
# filter matrix by ticks to match dimension of cemat and dyncc mat
def trimbytick(args):
    assert len(args) == 4, 'Usage: python utils_mat.py filterbytick matfile fulltickfile trimtickfile outfile'
    matfile = args[0]
    fulltickfile = args[1]
    trimtickfile = args[2]
    outfile = args[3]

    mat = np.loadtxt(matfile)

    fullticklist = cp.loadtuples(fulltickfile)[0]
    trimticklist = cp.loadtuples(trimtickfile)[0]
    assert len(fullticklist) > len(trimticklist), 'fulltick is smaller than trimtick.'

    outmat = _trimbytick(mat, fullticklist, trimticklist)
    np.savetxt(outfile, outmat)


def _fillbyticksep(mat, t, s, v):
    n = len(t)
    for i in range(n):
        for j in range(i,n):
            if t[j]-t[i] <= s:
                mat[i,j] = mat[j,i] = v
    return mat

# thrsholding matrix by tick separations
# for example filtering cemat by large sequential separations
def fillbyticksep(args):
    assert len(args) == 5, 'Usage:python utils_mat.py fillbyticksep inmatfile tickfile separation{int value} value outfile'

    inmatfile = args[0]
    tickfile = args[1]
    s = int(args[2])
    v = float(args[3])
    outfile = args[4]

    mat = np.loadtxt(inmatfile)
    t = np.loadtxt(tickfile).astype(int)
    n = len(t)
    assert mat.shape[0] == n, 'matrix dimension and number of ticks not match.'

    outmat = _fillbyticksep(mat, t, s, v)
    np.savetxt(outfile, outmat, fmt='%.3f')
    cp._info('save filtered matrix to %s' % outfile)


def foo(args):
    print(len(args))

if __name__=='__main__':
    cp.dispatch(__name__)