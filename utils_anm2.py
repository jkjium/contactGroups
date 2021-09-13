import numpy as np
import commp as cp

from protein import protein

class anm:
    def __init__(self, atoms, cutoff=12.0, gamma=1.0, ctmat=None):
        assert atoms!=None, 'No atoms info'
        self.atoms = atoms
        self.cutoff = cutoff

        # get contact matrix
        if ctmat==None:
            self.atoms = atoms
            self.coords = np.array([[at.x, at.y, at.z] for at in self.atoms])
            # squared pairwise distance matrix by numpy broadcasting
            self.spdmat = ((self.coords[:, :, None] - self.coords[:, :, None].T) ** 2).sum(1)
            # contact matrix
            self.ctmat = (self.spdmat <= self.cutoff**2).astype(int)
            np.fill_diagonal(self.ctmat, 0)
        else:
            self.ctmat = ctmat
            cp._info('External contact map used.')
        n = len(self.atoms)
        assert n == self.ctmat.shape[0], 'atoms and contact matrix have different length'

        self.hessmat = np.zeros((3*n,3*n),dtype = float)
        # filling hessian by iterating all contact pairs
        x,y =np.where(np.triu(self.ctmat)==1)
        for i,j in zip(x,y):
            d = self.coords[i] - self.coords[j]
            h3x3 = -gamma*np.outer(d,d)/self.spdmat[i,j]
            self.hessmat[3*i:3*i+3, 3*j:3*j+3] = self.hessmat[3*j:3*j+3, 3*i:3*i+3] = h3x3
            # accumulate diagonal
            self.hessmat[3*i:3*i+3, 3*i:3*i+3] -= h3x3
            self.hessmat[3*j:3*j+3, 3*j:3*j+3] -= h3x3

        # calculating modes (sorted by default)
        e, v = np.linalg.eigh(self.hessmat)
        self.modes = v
        self.frequencies = e
        
    def calcdyncc(self, mode_list):
        pass
    
    def calcdistflucts(self, mode_list):
        pass
    
    def calcmsf(self, mode_list):
        pass




def foo(args):
    p = protein('1ggg_a.pdb')
    a= anm(p.ca)
    #np.savetxt('t.hess.out', a.hessmat, fmt='%.3f')
    #np.savetxt('t.ctmat.out', a.ctmat, fmt='%d')
    print(a.frequencies.shape)
    for i in range(len(a.frequencies)):
        print('%d %.4f' % (i,a.frequencies[i]))

if __name__=='__main__':
    cp.dispatch(__name__)