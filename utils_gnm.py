import numpy as np
import commp as cp

from protein import protein
'''
calculate gnm freqieunces, modes, contact matrix, squared pdist matrix,
     cross-correlations, distance fluctuations, squared fluctuations
'''
class gnm:
    def __init__(self, pdbfile, ca=True, cutoff=10, gamma=1.0):
        self.p = protein(pdbfile)
        self.atoms = self.p.ca if ca == True else self.p.atoms
        self.cutoff = cutoff
        self.gamma = gamma
        self.coords = np.array([[at.x, at.y, at.z] for at in self.atoms])

        # squared pairwise distance matrix
        self.spdmat = ((self.coords[:, :, None] - self.coords[:, :, None].T) ** 2).sum(1)
        # contact matrix
        self.ctmat = (self.spdmat <= self.cutoff**2).astype(int)
        np.fill_diagonal(self.ctmat, 0)
        D = np.diag(np.sum(self.ctmat, axis=0))
        # laplacian
        self.kirchhoff = gamma * (D - self.ctmat)

        # modes
        e, v = np.linalg.eigh(self.kirchhoff)
        self.modes = v
        self.frequencies = e

    '''
    calc{funs}: mode_list must contain 1-based indices to exclude mode[0]
    '''
    # dynamic cross-correlation
    def calcdyncc(self, mode_list):
        v = self.modes[:,mode_list]
        d = np.diag(1/self.frequencies[mode_list])
        covariance = np.dot(v, np.dot(d, v.T))
        # normalization
        diag = np.power(covariance.diagonal(), 0.5)
        diagouter = np.outer(diag, diag)
        return cp.div0(covariance, diagouter)

    # return distance fluctuations matrix
    def calcdistflucts(self, mode_list):
        ccmat = self.calcdyncc(mode_list)
        cc_diag = np.diag(ccmat).reshape(-1,1)
        return cc_diag.T + cc_diag -2.*ccmat

    # squared fluctuations    
    def calcmsf(self, modelist):
        v = self.modes[:,modelist]
        diagvar = np.diag(1/self.frequencies[modelist])
        return np.dot(v * v, diagvar).sum(axis=1)

# output all ..
# calculate 20 modes by default
def analysis(args):
    assert len(args) >=4, 'Usage: python utils_gnm.py analysis pdbfile cutoff gamma is_ca (mode_list{1,2,3,...])'
    pdbfile = args[0]
    cutoff = float(args[1])
    gamma = float(args[2])
    is_ca = bool(args[3]) # or str in next version 
    mode_list = map(int, args[4].split(',')) if len(args) == 5 else list(range(1,20))

    g = gnm(pdbfile, is_ca, cutoff, gamma)
    dcmat = g.calcdyncc(mode_list)
    dfmat = g.calcdistflucts(mode_list)
    msf = g.calcmsf(mode_list)

    outadjmat = '%s.adjmat.txt' % pdbfile
    np.savetxt(outadjmat, g.ctmat, fmt='%d')

    outdcmat = '%s.dcmat.txt' % pdbfile
    np.savetxt(outdcmat, dcmat, fmt='%.3f')

    outdfmat = '%s.dfmat.txt' % pdbfile
    np.savetxt(outdfmat, dfmat, fmt='%.3f')

    outmsf = '%s.msf.txt' % pdbfile
    np.savetxt(outmsf, msf, fmt='%.3f')

    outmodes = '%s.modes.txt' % pdbfile
    np.savetxt(outmodes, g.modes, fmt='%.3f')

    outfreqs = '%s.freqs.txt' % pdbfile
    np.savetxt(outfreqs, g.frequencies, fmt='%.3f')

    cp._info('save %s.{adjmat, dcmat, dfmat, msf, modes, freqs}.txt' % pdbfile)
    

def foo(args):
    #g = gnm('t1.pdb')
    g = gnm('1aar_a.pdb')
    print(g.kirchhoff)
    mlist = [1,2,3]
    dyncc = g.calcdyncc(mlist)
    print('dyncc:')
    print(dyncc.round(3))
    print('distflucts:')
    print(g.calcdistflucts(mlist).round(3))
    print('---------------')
    print(g.calcmsf(mlist).round(3))

if __name__=='__main__':
    cp.dispatch(__name__)
