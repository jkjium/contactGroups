##############################################################
# acknowledge Yuan Wang as the author of anm implementation 
##############################################################

import commp as cp
import numpy as np

###################################
###### Simple Math Functions ######
###################################

def dis(cor_a, cor_b):
    '''
    returns distance between two points in 3 dimensions given the coordinates of the two points
    '''    
    summ = 0    
    for i in range (0,3):
        summ = summ + (cor_a[i]-cor_b[i]) ** 2
        
    distance = summ ** 0.5
    return distance

def mat2vec(np_matrix):
    '''
    transfers a matrix to a vector
    '''    
    vector = np.asarray(np_matrix).reshape(-1)
    return vector
    
def vec2mat(np_vector):
    '''
    transfers a vector to a 3D matrix
    '''
    
    length = len(np_vector)
    if length%3 != 0:
        raise ValueError('vector length cannot be divided by 3')
    else:
        matrix = np_vector.reshape((length/3,3))
    
    return matrix
    
##################################
########## ANM Class #############
##################################    

class ANM:
    def __init__(self, protein, cutoff = 12, useCA = True, method = 'cutoff', power = 4):
        '''
        initiation of class ANM
        '''        
        
        
        self.cx = None
        self.hess = None
        self.freq = None
        self.modes = None
        self.cutoff = cutoff
        self.e = None
        self.v = None
        self.proteinName = protein
        self.method = method
        self.power = power
        self.MSF = None
        self.arrows = list()
        global struc        
        try:
            struc = structure(protein)
        except IOError:
            if havePymol:
                obj_list = cmd.get_names()
                filename = protein + '_PyANM.pdb'
                if protein in obj_list:
                    cmd.save(filename, protein, state = -1 )
                else:
                    cmd.fetch(protein)
                    cmd.save(filename, protein, state = -1)
                struc = structure(filename)
                os.remove(filename)
            else:
                raise IOError('Failed to find object %s both in PDB and local directory' %protein)
              
        
        if useCA:
            struc.CA_only()
            if struc.isCA():
                pass
            else:
                raise IOError('Failed to extract all alpha carbons from %s' %protein)
        self.n_atoms = struc.length
        
        
    
    def getCX(self):
        if self.cx is None:
            cp._info("Contact Matrix not found, Rebuilding now...") 
            self.buildCX()
        
        return self.cx
            
    def resetCX(self):
        self.cx = None
        
    def getHess(self):
        if self.hess is None:
            cp._info("Hessian Matrix not found, Rebuilding now...")
            self.buildHess()
        
        return self.hess
        
    def resetHess(self):
        self.hess = None
        
    def getMSF(self):
        if self.MSF is None:
            self.calcMSF()
            
        return self.MSF
        
    def resetMSF(self):
        self.MSF = None
        
    def getE(self):
        if self.e is None:
            self.calcModes()
            
        return self.e
    
    def resetE(self):
        self.e = None
        
    def getV(self):
        if self.v is None:
            self.calcModes()
            
        return self.v
        
    def resetV(self):
        self.v = None
        
    def buildCX(self):
        '''
        build contact matrix
        '''
        if not isinstance(struc.cor, np.ndarray):        
            raise ValueError('Structure coordinates are not numpy arrays')
        if struc.cor.shape[1] != 3:
            raise ValueError('Dimension of coordinates does not equal 3')
            
        self.contactPair_x = list()
        self.contactPair_y = list()    
        length = struc.cor.shape[0]    
        self.cx = np.zeros((length,length),dtype = float)
        
        if self.method == 'cutoff':
            for i in range (1,length):
                for j in range (0,i):
                    temp_distance = dis(struc.cor[i,:], struc.cor[j,:])
                    if temp_distance < self.cutoff:
                        self.cx[i,j] = 1
                        self.contactPair_x.append(i)
                        self.contactPair_y.append(j)
        elif self.method == 'pf':
            for i in range (1,length):
                for j in range (0,i):
                    temp_distance = dis(struc.cor[i,:], struc.cor[j,:])
                    self.cx[i,j] = (1/temp_distance)**self.power
                    if temp_distance <= self.cutoff:
                        self.contactPair_x.append(i)
                        self.contactPair_y.append(j)
        else:
            raise IOError('Wrong ANM Method, should be either cutoff or pf')
                    
    def buildHess(self):
        '''
        build hessian matrix
        '''
        self.cx = self.getCX()
        length = self.cx.shape[0]
        contactNumber = len(self.contactPair_x)
        self.hess = np.zeros((3*length,3*length),dtype = float)
        for i in range (0,contactNumber):
            x = self.contactPair_x[i]
            y = self.contactPair_y[i]
            dR = struc.cor[y,:] - struc.cor[x,:]
            dR = dR/np.linalg.norm(dR)
            K33 = -dR[np.newaxis,:].T*dR
            self.hess[3*x:3*x+3,3*y:3*y+3] = K33
            self.hess[3*y:3*y+3,3*x:3*x+3] = K33
            self.hess[3*x:3*x+3,3*x:3*x+3] = self.hess[3*x:3*x+3,3*x:3*x+3] - K33
            self.hess[3*y:3*y+3,3*y:3*y+3] = self.hess[3*y:3*y+3,3*y:3*y+3] - K33
            
    def calcModes(self,numeig = None):
        '''
        calculat modes for ANM
        
        numeig: number of eigen values calculated
                default value is none
                all modes are calculated
        '''
        self.hess = self.getHess()
        if not isinstance(self.hess, np.ndarray):        
            raise ValueError('Hessian Matrix is not in the type of numpy arrays')
            
        if numeig is None:
            numeig = struc.cor.shape[0]*3
        '''
        try: 
            import scipy.sparse.linalg
            haveSP=True
        except ImportError:
            haveSP=False
        '''
        haveSP=False # scipy is slower and inaccurate
        if haveSP:
            try:
                [self.e,self.v] = scipy.sparse.linalg.eigsh(self.hess, numeig, which = 'SM', maxiter = 5000)
            except ValueError:
                [self.e,self.v] = scipy.sparse.linalg.eigsh(self.hess, numeig-1, which = 'SM', maxiter = 5000)
        else:
            [self.e,self.v] = np.linalg.eig(self.hess)
            self.v = self.v[:,self.e.argsort()]
            self.e.sort()
            self.v = self.v[:,0:numeig]
            self.e = self.e[0:numeig]
            
        if not self.freqChecker():
            cp._info("Not all eigenvalues converged, results might be inaccurate")
            print(self.e[0:7])
            
    
    def freqChecker(self):
        '''
        checks if the first 6 modes have a frequency of 0
        '''
        flag = False
        if sum(self.e[0:6]) < 1e-10 and sum(self.e[0:7]) > 1e-10:
            flag = True
        
        return flag
    
    def modeAnimator(self, modes = [1,2,3], useAllAtom = False, scaler = 5):
        '''
        make movies for different modes
        
        modes:        mode numbers to make movie for
                      default is first 3 modes
        scaler:       motion scales for movie
                      default is 5
        '''
        self.v = self.getV()
        self.movieFiles = {}
        numberTracker = 0
        if not useAllAtom:
            for i in modes:
                movieFile = self.proteinName.split('.')[0] + '_mode' + str(i) + '.pdb'
                cp._info('save to movie file: %s' % movieFile)
                self.movieFiles[numberTracker] = movieFile
                numberTracker += 1
                outfile = open(movieFile,'w')
                
                for j in range(11):
                    tempCor = struc.cor + scaler*j*vec2mat(self.v[:,i+5]) # +5 skip the first 6 modes
                    outfile.write('MODEL  %4d\n' %(j+1))
                    for k in range(struc.length):
                        outfile.write('%-6s%5d  %-4s%3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f           %-4s\n' %(struc.line_type[k], struc.atom_number[k], struc.atom_name[k], struc.residue_name[k], struc.chain_id[k], struc.residue_number[k],tempCor[k,0], tempCor[k,1], tempCor[k,2], struc.occupancy[k], struc.b_factor[k], struc.element[k]))
                    outfile.write('ENDMDL\n')
                
                for j in range(10):
                    tempCor = struc.cor + scaler*(10-j)*vec2mat(self.v[:,i+5])
                    outfile.write('MODEL  %4d\n' %(j+12))
                    for k in range(struc.length):
                        outfile.write('%-6s%5d  %-4s%3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f           %-4s\n' %(struc.line_type[k], struc.atom_number[k], struc.atom_name[k], struc.residue_name[k], struc.chain_id[k], struc.residue_number[k],tempCor[k,0], tempCor[k,1], tempCor[k,2], struc.occupancy[k], struc.b_factor[k], struc.element[k]))
                    outfile.write('ENDMDL\n')
                
                for j in range(11):
                    tempCor = struc.cor - scaler*j*vec2mat(self.v[:,i+5])
                    outfile.write('MODEL  %4d\n' %(j+22))
                    for k in range(struc.length):
                        outfile.write('%-6s%5d  %-4s%3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f           %-4s\n' %(struc.line_type[k], struc.atom_number[k], struc.atom_name[k], struc.residue_name[k], struc.chain_id[k], struc.residue_number[k],tempCor[k,0], tempCor[k,1], tempCor[k,2], struc.occupancy[k], struc.b_factor[k], struc.element[k]))
                    outfile.write('ENDMDL\n')
                    
                for j in range(10):
                    tempCor = struc.cor - scaler*(10-j)*vec2mat(self.v[:,i+5])
                    outfile.write('MODEL  %4d\n' %(j+33))
                    for k in range(struc.length):
                        outfile.write('%-6s%5d  %-4s%3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f           %-4s\n' %(struc.line_type[k], struc.atom_number[k], struc.atom_name[k], struc.residue_name[k], struc.chain_id[k], struc.residue_number[k],tempCor[k,0], tempCor[k,1], tempCor[k,2], struc.occupancy[k], struc.b_factor[k], struc.element[k]))
                    outfile.write('ENDMDL\n')
                
                outfile.close()
                
    def calcMSF(self):
        '''
        calculate the Mean Square Fluctuations from ANM
        '''
        self.cx = self.getCX()
        self.hess = self.getHess()
        
        if self.e is None:
            self.calcModes()
            
        self.MSF = np.zeros((1,struc.length), dtype = float)
        
        cv1 = range(0,struc.length*3,3)
        cv2 = range(1,struc.length*3,3)
        cv3 = range(2,struc.length*3,3)
        
        for i in range(6,struc.length*3):
            v_temp = self.v[:,i]**2
            self.MSF = self.MSF + (1/(self.e[i])) * (v_temp[cv1]+v_temp[cv2]+v_temp[cv3])
        
        T  = 300
        Kb = 1.38065e-23
        Gamma = 1e-20
        pi = 3.1415926
         
        self.MSF =  (8/3)*(pi**2)*Kb*(T/Gamma)*self.MSF
        
    def arrowGenerator(self, mode = 1, color = [1,1,1], size = 5):
        '''
        generate arrow objects for pymol
        
        mode:   the number of mode to draw arrows for
                default is the first mode
        color:  rgb color for the arrows
                default is white
        size:   size for the arrows
                default is 5
        '''
        self.v = self.getV()
        self.arrows = list()
        
        for i in range(struc.length):
            tail = [CYLINDER, struc.cor[i,0], struc.cor[i,1], struc.cor[i,2],\
                    struc.cor[i,0]+size*16*self.v[3*i,mode+5], struc.cor[i,1]+size*16*self.v[3*i+1,mode+5], struc.cor[i,2]+size*16*self.v[3*i+2,mode+5],\
                    0.3, color[0], color[1], color[2], color[0], color[1], color[2]]
            self.arrows.extend(tail)
            
            head = [CONE, struc.cor[i,0]+size*16*self.v[3*i,mode+5], struc.cor[i,1]+size*16*self.v[3*i+1,mode+5], struc.cor[i,2]+size*16*self.v[3*i+2,mode+5],\
                   struc.cor[i,0]+size*20*self.v[3*i,mode+5], struc.cor[i,1]+size*20*self.v[3*i+1,mode+5], struc.cor[i,2]+size*20*self.v[3*i+2,mode+5],\
                   1, 0.0, color[0], color[1], color[2], color[0], color[1], color[2], 1.0, 1.0]
            self.arrows.extend(head)
        #cmd.load_cgo(self.arrows,'test')
        
##################################
###### PDB Parsing Class #########
##################################

class structure:
    def __init__(self, filename, includeHET = False):
        '''
        A pdb parser for formal pdb format
        
        filename:   the pdb file to be parsed
        includeHET: if HETATM in pdb should be included, default is false (do not include)
        '''
        pdbfile = open(filename,'r')
        self.line_type = list()
        self.atom_number = list()
        self.atom_name = list()
        self.residue_name = list()
        self.chain_id = list()
        self.residue_number = list()
        self.cor_x = list()
        self.cor_y = list()
        self.cor_z = list()
        self.occupancy = list()       
        self.b_factor = list()
        self.element = list()
        
        for line in pdbfile:
            if includeHET:
                if (line[0:4] == 'ATOM') | (line[0:6] == 'HETATM'):
                    self.line_type.append(line[0:6].strip())                
                    self.atom_number.append(int(line[7:11].strip()))
                    self.atom_name.append(line[12:16].strip())
                    self.residue_name.append(line[17:21].strip())
                    self.chain_id.append(line[21])
                    self.residue_number.append(int(line[22:26].strip()))
                    self.cor_x.append(float(line[30:38].strip()))
                    self.cor_y.append(float(line[38:46].strip()))
                    self.cor_z.append(float(line[46:54].strip()))
                    self.occupancy.append(float(line[54:60].strip()))
                    self.b_factor.append(float(line[60:66].strip()))
                    self.element.append(line[77:80].strip())
            else:
                if (line[0:4] == 'ATOM'):
                    self.line_type.append(line[0:6].strip())                
                    self.atom_number.append(int(line[7:11].strip()))
                    self.atom_name.append(line[12:16].strip())
                    self.residue_name.append(line[17:21].strip())
                    self.chain_id.append(line[21])
                    self.residue_number.append(int(line[22:26].strip()))
                    self.cor_x.append(float(line[30:38].strip()))
                    self.cor_y.append(float(line[38:46].strip()))
                    self.cor_z.append(float(line[46:54].strip()))
                    self.occupancy.append(float(line[54:60].strip()))
                    self.b_factor.append(float(line[60:66].strip()))
                    self.element.append(line[77:80].strip())
        
        self.length = len(self.atom_name)
        self.cor = np.zeros((self.length,3),dtype=float)
        self.cor[:,0] = self.cor_x;
        self.cor[:,1] = self.cor_y;
        self.cor[:,2] = self.cor_z;
        pdbfile.close()
    
    def CA_only(self):
        '''
        parse only the alpha carbons
        '''
        temp_line_type = list()        
        temp_atom_number = list()
        temp_atom_name = list()
        temp_residue_name = list()
        temp_chain_id = list()
        temp_residue_number = list()
        temp_occupancy = list()
        temp_b_factor = list()
        temp_element = list()
        
        ca_index = list()
        for i in range(self.length):
            if self.atom_name[i] == 'CA':
                ca_index.append(i)
        
        temp_cor = self.cor[ca_index,:]
        ca_index.reverse()
        while ca_index:
            temp = ca_index.pop()
            temp_line_type.append(self.line_type[temp])            
            temp_atom_number.append(self.atom_number[temp])
            temp_atom_name.append(self.atom_name[temp])
            temp_residue_name.append(self.residue_name[temp])
            temp_chain_id.append(self.chain_id[temp])
            temp_residue_number.append(self.residue_number[temp])
            temp_occupancy.append(self.occupancy[temp])
            temp_b_factor.append(self.b_factor[temp])
            temp_element.append(self.element[temp])
        
        self.line_type = temp_line_type        
        self.atom_number = temp_atom_number
        self.atom_name = temp_atom_name
        self.residue_name = temp_residue_name
        self.chain_id = temp_chain_id
        self.residue_number = temp_residue_number
        self.occupancy = temp_occupancy
        self.b_factor = temp_b_factor
        self.element = temp_element
        self.cor = temp_cor
        self.length = len(self.atom_name)
        
    def isCA(self):
        '''
        returns true if the structure contains only alpha carbons
        '''
        ca_flag = True        
        for i in range(self.length):
            if self.atom_name[i] == 'CA':
                continue
            else:
                ca_flag = False
                break
        return ca_flag
        
    def writePDB(self,outname):
        '''
        write structure back to pdb format
        
        outname: output filename
        '''
        outfile = open(outname,'w')
        for i in range(self.length):
           outfile.write('%-6s%5d  %-4s%3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f           %-4s\n' %(self.line_type[i], self.atom_number[i], self.atom_name[i], self.residue_name[i], self.chain_id[i], self.residue_number[i],self.cor[i,0], self.cor[i,1], self.cor[i,2], self.occupancy[i], self.b_factor[i], self.element[i]))
       
        outfile.close()


####################################
# executable procedures
####################################

def _ANMbuilder(pdbfile, cutoff):
    anm = ANM(pdbfile, cutoff)
    anm.buildHess()
    anm.calcModes()
    if anm.freqChecker():
        cp._info('ANM Built successfully')
    else:
        cp._info('Warning, Not all of first 6 eigenvalues are zero')
    return anm


def _inverseHessian(modes, eigs):
    n_modes = modes.shape[1] if len(modes.shape) > 1 else 1
    dim = len(modes)
    print n_modes, dim
    outersum = np.zeros((dim, dim))
    for i in range(n_modes):
        outersum+=(1/eigs[i])*np.outer(modes[:,i],modes[:,i])
    return outersum


# calculate cross-correlations from given modes and eigenvalues    
# cc_ij = tr(H^-1_ij) / sqrt(tr(H^-1_ii) * tr(H^-1_jj))
#    ccmat = _crosscorrelation(anm.v[:,mode_correlation_list+5], anm.e[mode_correlation_list+5])
def _crosscorrelation_formula(modes, eigs):
    assert len(modes) % 3 ==0, 'Error: mode length cannot be divided by 3'
    iH = _inverseHessian(modes, eigs)
    dim = len(modes)/3
    ccmat = np.zeros((dim, dim))
    for i in range(dim):
        for j in range(dim):
            trij = np.trace(iH[i*3:i*3+3, j*3:j*3+3])
            trii = np.trace(iH[i*3:i*3+3, i*3:i*3+3])
            trjj = np.trace(iH[j*3:j*3+3, j*3:j*3+3])
            ccmat[i,j] = trij / np.sqrt(trii*trjj)
    return ccmat


def _div0(a, b):
    """ Performs ``true_divide`` but ignores the error when division by zero 
    (result is set to zero instead). """

    from numpy import errstate, true_divide, isfinite, isscalar
    
    with errstate(divide='ignore', invalid='ignore'):
        c = true_divide(a, b)
        if isscalar(c):
            if not isfinite(c):
                c = 0
        else:
            c[~isfinite(c)] = 0.  # -inf inf NaN
    return c

def _par_crosscorrelations(queue, n_atoms, array, variances, indices):
    """Calculate covariance-matrix for a subset of modes."""

    n_modes = len(indices)
    arvar = (array[:, indices] * variances[indices]).T.reshape((n_modes,
                                                                n_atoms, 3))
    array = array[:, indices].T.reshape((n_modes, n_atoms, 3))
    covariance = np.tensordot(array.transpose(2, 0, 1),
                              arvar.transpose(0, 2, 1),
                              axes=([0, 1], [1, 0]))
    queue.put(covariance)

# performance version of cc calculation 
def _crosscorrelation(anm, mode_list, n_cpu=1):
    n_modes = len(mode_list)
    n_atoms = anm.n_atoms
    indices = mode_list

    array = anm.v.copy()
    variances = 1.0/anm.e # variance is the inverse of the eigenvalue

    if n_cpu == 1:
        s = (n_modes, n_atoms, 3)
        arvar = (array[:, indices]*variances[indices]).T.reshape(s)
        array = array[:, indices].T.reshape(s)
        covariance = np.tensordot(array.transpose(2, 0, 1),
                                  arvar.transpose(0, 2, 1),
                                  axes=([0, 1], [1, 0]))
    else:
        import multiprocessing
        n_cpu = min(multiprocessing.cpu_count(), n_cpu)
        queue = multiprocessing.Queue()
        size = n_modes / n_cpu
        for i in range(n_cpu):
            if n_cpu - i == 1:
                indices = modes.indices[i*size:]
            else:
                indices = modes.indices[i*size:(i+1)*size]
            process = multiprocessing.Process(
                target=_par_crossCorrelations,
                args=(queue, n_atoms, array, variances, indices))
            process.start()
        while queue.qsize() < n_cpu:
            time.sleep(0.05)
        covariance = queue.get()
        while queue.qsize() > 0:
            covariance += queue.get()

    # if norm:
    diag = np.power(covariance.diagonal(), 0.5)
    D = np.outer(diag, diag)
    covariance = _div0(covariance, D)

    return covariance


# analysis outputs: 
#   _frequencies.txt: eigenvalues (single column file, len = n_atoms*3)
#   _modes.txt: eigenvectors (each column is a mode, n_atoms*3 x n_atoms*3)
#   _contactmatrix.txt: a square matrix (n_atoms x n_atoms)
#   _Hessian.txt: a square matrix (n_atoms*3 x n_atoms*3)
#   _MSF.txt: single row file  (len = n_atoms)
#   _ccmat: cross-correlation matrix (n_atoms x n_atoms)
def anmanalysis(args):
    assert len(args) == 5, 'Usage: python anm pdbfile cutoff modes_for_animation {"1,2,3"} animate_scaler modes_for_correlation {"1,2,3"}'
    pdbfile = args[0]
    cutoff = float(args[1])
    mode_animate_list = list(map(int, args[2].split(','))) 
    animate_scaler = float(args[3])
    mode_correlation_list = np.array(list(map(int, args[4].split(',')))) # must be a list, even only one mode selected

    anm = _ANMbuilder(pdbfile, cutoff)

    np.savetxt(anm.proteinName + '_frequencies.txt', anm.e, fmt = '%.3f')
    np.savetxt(anm.proteinName + '_modes.txt', anm.v, fmt = '%.3f')
    np.savetxt(anm.proteinName + '_contactMatrix.txt', anm.getCX(), fmt = '%d')
    #np.savetxt(anm.proteinName + '_contactMatrix.txt', anm.getCX(), fmt = '%.3e') # for method='pf'
    np.savetxt(anm.proteinName + '_Hessian.txt', anm.getHess(), fmt = '%.3f')
    np.savetxt(anm.proteinName + '_MSF.txt', anm.getMSF(), fmt = '%.3f')
    cp._info('Data (frequencies, modes, contactMatrix, Hessian, MSF) saved.')

    anm.modeAnimator(mode_animate_list, scaler=animate_scaler)
    cp._info('%d animation movies saved.' % len(mode_animate_list))
    
    #ccmat = _crosscorrelation(anm.v[:,mode_correlation_list+5], anm.e[mode_correlation_list+5])
    ccmat = _crosscorrelation(anm, mode_correlation_list+5)
    print ccmat.round(3)
    np.savetxt(anm.proteinName + '_ccmat.txt', ccmat, fmt = '%.3f')
    cp._info('Cross-Correlation matrix of %d modes : %s_ccmat.txt saved.' % (len(mode_correlation_list), anm.proteinName))


if __name__ == '__main__':
    cp.dispatch(__name__)