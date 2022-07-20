import numpy as np
import commp as cp

from protein import protein
'''
calculate gnm freqieunces, modes, contact matrix, squared pdist matrix,
     cross-correlations, distance fluctuations, squared fluctuations
'''
class gnm:
    def __init__(self, atoms=None, cutoff=10, gamma=1.0, ctmat=None):
        self.cutoff = cutoff
        self.gamma = gamma

        if ctmat is None:
            assert atoms!=None
            self.atoms = atoms
            self.coords = np.array([[at.x, at.y, at.z] for at in self.atoms])
            # squared pairwise distance matrix
            self.spdmat = ((self.coords[:, :, None] - self.coords[:, :, None].T) ** 2).sum(1)
            # contact matrix
            self.ctmat = (self.spdmat <= self.cutoff**2).astype(int)
        else:
            self.ctmat = ctmat
            #cp._info('External contact map used.')

        np.fill_diagonal(self.ctmat, 0)
        D = np.diag(np.sum(self.ctmat, axis=0))
        # laplacian
        self.kirchhoff = gamma * (D - self.ctmat)

        # modes
        e, v = np.linalg.eigh(self.kirchhoff)
        self.modes = v # when using skip the first one
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

    p = protein(pdbfile)
    g = gnm(atoms=p.ca, cutoff=cutoff)
    outticks = '%s.gnm.tick' % pdbfile
    with open(outticks, 'w') as fout:
        fout.write('%s\n' % (' '.join([str(a.resSeq) for a in p.ca])))

    dcmat = g.calcdyncc(mode_list)
    dfmat = g.calcdistflucts(mode_list)
    msf = g.calcmsf(mode_list)

    outadjmat = '%s.gnm.adjmat.txt' % pdbfile
    np.savetxt(outadjmat, g.ctmat, fmt='%d')

    outdcmat = '%s.gnm.dcmat.txt' % pdbfile
    np.savetxt(outdcmat, dcmat, fmt='%.3f')

    outdfmat = '%s.gnm.dfmat.txt' % pdbfile
    np.savetxt(outdfmat, dfmat, fmt='%.3f')

    outmsf = '%s.gnm.msf.txt' % pdbfile
    np.savetxt(outmsf, msf, fmt='%.3f')

    outmodes = '%s.gnm.modes.txt' % pdbfile
    np.savetxt(outmodes, g.modes, fmt='%.8f')

    outfreqs = '%s.gnm.freqs.txt' % pdbfile
    np.savetxt(outfreqs, g.frequencies, fmt='%.8f')

    cp._info('save %s.gnm.{tick, adjmat, dcmat, dfmat, msf, modes, freqs}.txt' % pdbfile)


# gnm using external contact map
def analysis_ctmat(args):
    assert len(args) ==1, 'Usage: python utils_gnm.py analysis ctmatfile outprefix'
    ctmatfile = args[0]

    #p = protein(pdbfile)
    ctmat = np.loadtxt(ctmatfile)
    g = gnm(ctmat=ctmat)
    '''
    dcmat = g.calcdyncc(mode_list)
    dfmat = g.calcdistflucts(mode_list)
    msf = g.calcmsf(mode_list)

    outadjmat = '%s.gnm.adjmat.txt' % pdbfile
    np.savetxt(outadjmat, g.ctmat, fmt='%d')
    outdcmat = '%s.gnm.dcmat.txt' % pdbfile
    np.savetxt(outdcmat, dcmat, fmt='%.3f')

    outdfmat = '%s.gnm.dfmat.txt' % pdbfile
    np.savetxt(outdfmat, dfmat, fmt='%.3f')

    outmsf = '%s.gnm.msf.txt' % pdbfile
    np.savetxt(outmsf, msf, fmt='%.3f')
    '''

    outmodes = '%s.gnm.modes.txt' % ctmatfile
    np.savetxt(outmodes, g.modes, fmt='%.8f')

    outfreqs = '%s.gnm.freqs.txt' % ctmatfile
    np.savetxt(outfreqs, g.frequencies, fmt='%.8f')
    # ticks are in 11AS_B-PF03590.dcaper.ccmat.tick from utils_ce.py cflat2ccmat
    # gnm tick needs to be converted with awk '{for(i=1;i<=NF;i++) print $i}' 11AS_B-PF03590.pdb.gnm.tick > new.tick
    cp._info('save %s.gnm.{adjmat, dcmat, dfmat, msf, modes, freqs}.txt' % ctmatfile)

# gnm using correlation matrix with ctmat cutoff
def analysis_cemat0(args):
    assert len(args) ==2, 'Usage: python utils_gnm.py analysis ctmatfile outprefix'
    cematfile = args[0]
    cutoff = float(args[1])
    outprefix = '%s.c%.4f' % (cematfile, cutoff)

    #p = protein(pdbfile)
    cemat = np.loadtxt(cematfile)
    ctmat = (cemat <= cutoff).astype(int)

    g = gnm(ctmat=ctmat)

    outadjmat = '%s.gnm.adjmat.txt' % outprefix
    np.savetxt(outadjmat, g.ctmat, fmt='%d')
    '''
    dcmat = g.calcdyncc(mode_list)
    dfmat = g.calcdistflucts(mode_list)
    msf = g.calcmsf(mode_list)

    outdcmat = '%s.gnm.dcmat.txt' % pdbfile
    np.savetxt(outdcmat, dcmat, fmt='%.3f')

    outdfmat = '%s.gnm.dfmat.txt' % pdbfile
    np.savetxt(outdfmat, dfmat, fmt='%.3f')

    outmsf = '%s.gnm.msf.txt' % pdbfile
    np.savetxt(outmsf, msf, fmt='%.3f')
    '''

    outmodes = '%s.gnm.modes.txt' % outprefix
    np.savetxt(outmodes, g.modes, fmt='%.8f')

    outfreqs = '%s.gnm.freqs.txt' % outprefix
    np.savetxt(outfreqs, g.frequencies, fmt='%.8f')
    # ticks are in 11AS_B-PF03590.dcaper.ccmat.tick from utils_ce.py cflat2ccmat
    # gnm tick needs to be converted with awk '{for(i=1;i<=NF;i++) print $i}' 11AS_B-PF03590.pdb.gnm.tick > new.tick
    cp._info('save %s.gnm.{adjmat, modes, freqs}.txt' % outprefix)


# gnm using correlation matrix with ctmat cutoff
# cutoff_start: percentile cutoff {0.04}. Iterate start from this value until a connected graph is found
# offset: construct ctmat using the minimum cutoff_start + offset that makes the new ctmat a connected graph
def analysis_cemat(args):
    assert len(args) ==4, 'Usage: python utils_gnm.py analysis cematfile cutoff_start offset searching_interval outprefix'
    cematfile = args[0]
    cutoff_start = float(args[1])
    offset = float(args[2])
    interval = float(args[3])

    cemat = np.loadtxt(cematfile)
    cutoff = 99
    # iterate until a connected graph is found
    #for c in range(cutoff_start, interval, 0.16):
    for c in np.arange(cutoff_start, 0.20, interval):
        print('iterating cutoff: %.2f' % c)
        ctmat = (cemat <= c).astype(int)
        g = gnm(ctmat=ctmat)
	#print(g.frequencies[0:5])
        if np.round(g.frequencies[1],8)!=0.0:
            cutoff = c
            break

    # apply ctmat cutoff offset
    if offset!=0:
        cutoff+=offset
        ctmat = (cemat <= (cutoff)).astype(int)
        g = gnm(ctmat=ctmat)


    outprefix = '%s.%.4f' % (cematfile, cutoff)
    # outputs
    # ctmat
    outadjmat = '%s.gnm.adjmat.txt' % (outprefix, cutoff)
    np.savetxt(outadjmat, g.ctmat, fmt='%d')

    # full modes matrix (including the first one)
    outmodes = '%s.gnm.modes.txt' % outprefix
    np.savetxt(outmodes, g.modes, fmt='%.8f')

    # frequencies (eigenvalues)
    outfreqs = '%s.gnm.freqs.txt' % outprefix
    np.savetxt(outfreqs, g.frequencies, fmt='%.8f')
    # ticks are in 11AS_B-PF03590.dcaper.ccmat.tick from utils_ce.py cflat2ccmat
    # gnm tick needs to be converted with awk '{for(i=1;i<=NF;i++) print $i}' 11AS_B-PF03590.pdb.gnm.tick > new.tick
    cp._info('save %s.gnm.{adjmat, modes, freqs}.txt' % outprefix)


    

def foo(args):
    p = protein('1ggg_a.pdb')
    g = gnm(atoms=p.ca)
    for i in range(len(g.frequencies)):
        print('%d %.4f' %  (i, g.frequencies[i]))
    '''
    print(g.kirchhoff)
    mlist = [1,2,3]
    dyncc = g.calcdyncc(mlist)
    print('dyncc:')
    print(dyncc.round(3))
    print('distflucts:')
    print(g.calcdistflucts(mlist).round(3))
    print('---------------')
    print(g.calcmsf(mlist).round(3))
    '''

if __name__=='__main__':
    cp.dispatch(__name__)
