import numpy as np
import commp as cp

# convert square matrix to flat file
def mat2flat(args):
    assert len(args) >= 3, 'Usage: python utils_mat.py mat2flat mat.file name.file outfile'
    matfile = args[0]
    namefile = args[1]
    outfile = args[2] 
    dc = args[3] if len(args) > 3 else ' '

    mat = np.loadtxt(matfile, delimiter=dc)
    n = mat.shape[0]
    name = cp.loadlines(namefile)

    outlist = []
    for i in range(n):
        for j in range(i+1,n):
            outlist.append('%s %s %.4f' % (name[i], name[j], mat[i,j]))

    with open(outfile, 'w') as fout:
        fout.write('%s\n' % ('\n'.join(outlist)))
    cp._info('save flat to %s' % outfile)
    

# calculate (accumulative) overlaps between two sets of vectors
# kuse: a vector of values for the number of modes to use in accumulative overlap 
def _modeoverlap(m1, m2, kuse):
    print(m1.shape)
    print(m2.shape)
    print(kuse)


# output mode (accumulative) overlaps
def modeoverlap(args):
    assert len(args) == 6, 'Usage: python utils_mat.py modeoverlap modes1.txt col_range1{1:10} mat2 col_range2{6:16} co_grid outprefix'
    matfile1 = args[0]
    mode_range1 = map(int, args[1].split(':'))
    matfile2 = args[2]
    mode_range2 = map(int, args[3].split(':'))
    cogrid = map(int, args[4].split(','))
    outprefix = args[5]

    mat1 = np.loadtxt(matfile1)
    mat2 = np.loadtxt(matfile1)

    # get focused modes to start comparison
    modes1 = mat1[:,mode_range1[0]:mode_range1[1]]
    modes2 = mat2[:,mode_range2[0]:mode_range2[1]]

    o, co, rmsip = _modeoverlap(modes1, modes2, cogrid)


# rmsip Amadei, A. et al. (1999)
def _rmsip(m1, m2):
    #print(m1.ndim) # single vector
    if m1.ndim==1: # single vector
	return np.dot(m1.T,m2)

    n = m1.shape[1]
    m1 *= 1 / (m1 ** 2).sum(0) ** 0.5
    m2 *= 1 / (m2 ** 2).sum(0) ** 0.5
    overlaps = np.dot(m1.T, m2)
    return np.sqrt(np.power(overlaps, 2).sum() / n)

# calculate rmsip
def rmsip(args):
    assert len(args) == 4, 'Usage: python utils_mat.py rmsip mode1.txt mode2.txt mode_range{1:10} outfile'
    matfile1 = args[0]
    matfile2 = args[1]
    m_range = args[2]
    mode_range = map(int, m_range.split(':'))
    outfile = args[3]

    # load mode matrices
    mat1 = np.loadtxt(matfile1)
    mat2 = np.loadtxt(matfile2)

    modes1 = mat1[:, mode_range[0]:mode_range[1]+1] if mode_range[0]< mode_range[1] else mat1[:,mode_range[0]]
    modes2 = mat2[:, mode_range[0]:mode_range[1]+1] if mode_range[0]< mode_range[1] else mat2[:,mode_range[0]]

    rmsip = _rmsip(modes1, modes2)
    outstr = 'rmsip %s %s %s %.4f' % (matfile1, matfile2, m_range, rmsip)
    print(outstr)
    with open(outfile, 'w') as fout:
        fout.write('%s\n' % outstr)

# calculate rmsip for a set of mode ranges
def rmsip_range(args):
    assert len(args) == 4, 'Usage: python utils_mat.py rmsip mode1.txt mode2.txt mode_range_set{1,2,3,5,10,20} outfile'
    matfile1 = args[0]
    matfile2 = args[1]
    mode_range = map(int, args[2].split(',')) # [1,2,3,5,10,20]
    outfile = args[3]

    # load mode matrices
    mat1 = np.loadtxt(matfile1)
    mat2 = np.loadtxt(matfile2)
    rmsiplist = []
    for r in mode_range:
        range_start = 1
        modes1 = mat1[:, range_start:r+1] if range_start < r else mat1[:,range_start]
        modes2 = mat2[:, range_start:r+1] if range_start < r else mat2[:,range_start]
        rmsiplist.append(_rmsip(modes1, modes2))
    outstr = '%s %s %s' % (matfile1, matfile2, ' '.join(['%.4f' % c for c in rmsiplist]))
    print(outstr)
    with open(outfile, 'w') as fout:
        fout.write('%s\n' % outstr)

# use trimtick to reduce fulltick matrix
# deal with dimension matching between dynmat and cemat
def _trimbytick(mat, fulltick, trimtick):
    trimidlist = [fulltick.index(i) for i in (set(fulltick) - set(trimtick))]
    outmat = np.delete(mat, trimidlist, axis=0) # trim rows
    outmat = np.delete(outmat, trimidlist, axis=1) # trim columns
    return outmat


# resi list from cemat is less than dyncc mat
# filter matrix by ticks to match dimension of cemat and dyncc mat
def trimbytick(args):
    assert len(args) == 4, 'Usage: python utils_mat.py trimbytick matfile fulltickfile trimtickfile outfile'
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

# thresholding matrix by tick separations
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
