import numpy as np
import commp as cp
from protein import protein
from atom import atom
from utils_pfammsa import pfammsa

# convert square matrix to vec3 flat
# for enm distance fluctuations matrix, dynamic cross-correlations matrix (should be symmetrical)
# input: pdbfile and matrix 
# output: upper diagonal (including diagonal) elements; a .vec3 flat file, r1 r2 value
def mat2resiflat(args):
    assert len(args) == 3, 'Usage: python proc_minip.py mat2resiflat pdbfile matfile outfile'
    pdbfile = args[0]
    matfile = args[1]
    outfile = args[2]

    p = protein(pdbfile)
    names = ['%d' % at.resSeq for at in p.ca]

    mat = np.loadtxt(matfile, delimiter=' ')
    assert len(names)==mat.shape[0], 'mismatched dimension: pdb: %d - mat: %d' % (len(names), mat.shape[0])

    outlist = ['%s %s %.3f' % (v[0], v[1], v[2]) for v in cp.mat2flat(mat,names)]
    with open(outfile, 'w') as fout:
        fout.write('%s\n' % '\n'.join(outlist))
    
    cp._info('save vec3 flat to %s' % outfile)

# print dot product and pearson correlations between $3 and $4
# used by correlations between dfmat and dca
def outcorrelation(args):
    assert len(args) == 2, 'Usage: python proc_minip.py printcorrelation *.dca.dfmat.vec3 outfile'
    d = np.array([[float(v[2]),float(v[3])] for v in cp.loadtuples(args[0])])
    nd0 = d[:,0]/np.linalg.norm(d[:,0])
    nd1 = d[:,1]/np.linalg.norm(d[:,1])
    pc = np.corrcoef(d[:,0], d[:,1])
    with open(args[1], 'w') as fout:
        fout.write('%s %.3f %.3f\n' % (args[0], np.dot(nd0,nd1), pc[0,1]))

# normalize ce and df scores in *.{dca}.dfmat.vec3 file
# (d-d.min)/d.max
def vec3norm(args):
    assert len(args) == 2, 'Usage: python proc_minip.py vec3norm 101m_A.pdb.PF00042.mip.dfmat.vec3 101m_A.pdb.PF00042.mip.dfmat.n.vec3'
    infile = args[0]
    outfile = args[1]
    resipairs = []
    scores = []
    for t in cp.loadtuples(infile):
        resipairs.append('%s %s' % (t[0], t[1])) # r1 r2
        scores.append([float(t[2]), float(t[3])])
    scoresnp = np.array(scores)
    # take min to remove negative values 
    #nce = 1.0*(scoresnp[:,0]-np.min(scoresnp[:,0]))/(np.max(scoresnp[:,0])-np.min(scoresnp[:,0]))
    #ndyn = 1.0*(scoresnp[:,1]-np.min(scoresnp[:,1]))/(np.max(scoresnp[:,1])-np.min(scoresnp[:,1]))
    nce = cp.normminmax(scoresnp[:,0])
    ndyn = cp.normminmax(scoresnp[:,1])

    outlist = ['%s %.8f %.8f' % (resipairs[i], nce[i], ndyn[i]) for i in range(len(resipairs))]
    with open(outfile, 'w') as fout:
        fout.write('%s\n' % '\n'.join(outlist))
    cp._info('save normalized vec3 to %s' % outfile)



# calculate weighed value density of 210 aa pairs from MSA {dca, mi, mip}
# output ce density and dfmat density
def aapaircorrdensity(args):
    assert len(args) == 5, 'Usage: python utils_minip.py *.{dca,mi,mip}.dfmat.nvec3 PF00000_p90.txt cflatfile weightfile outfile'
    # get 210 aa pairs (weighted) frequenccies
    dynfile = args[0]
    msafile = args[1]
    cflatfile = args[2]
    weightfile = args[3]
    outprefix = args[4]

    # load resipair2msaipair map from cflat
    # map[$2,$3] = $12,$13
    pmap = dict(('%s %s' % (ce[2],ce[3]), '%s %s' % (ce[12], ce[13])) for ce in cp.loadtuples(cflatfile))

    # load ce.dyn file
    # r1  r2  dca      distance_fluctuations
    # 593 594 0.342934 0.014
    corrdict = dict(('%s %s' % (c[0], c[1]), (float(c[2]), float(c[3]))) for c in cp.loadtuples(dynfile))

    # load msa
    pfm = pfammsa(msafile)
    w = np.loadtxt(weightfile)

    # init
    aa = cp.aat01
    aaidx =  ['%s %s' % (aa[i], aa[j]) for i in range(len(aa)) for j in range(len(aa))]
    freqdensity = np.zeros(400)
    cedensity = np.zeros(400)
    dyndensity = np.zeros(400)

    # loop through all residue pairs
    count=0
    for k in corrdict:
        # get msai pair
        plist = map(int, pmap[k].split(' '))
        # calcualte weighed 400 aa frequencies for current pairs of positions
        pfreqwdict = pfm.aapairfreqw(plist, w, aa)
        #print(pfreqwdict)

        # convert to np arrays to multiply ce score and dyn cc score
        pfreqwnp_o = np.array([pfreqwdict[aapair] for aapair in aaidx])
        #pfreqwnp = cp.normminmax(pfreqwnp_o)
        pfreqwnp = pfreqwnp_o

        '''
        # save each pair to a file for debug purpose
        outname = 't.s.%d.%d.%d.dcamat' % (count, plist[0], plist[1])
        outdca = pfreqwnp * corrdict[k][0]
        np.savetxt(outname, outdca.reshape(20,20), fmt='%.4f')

        outname = 't.s.%d.%d.%d.freqmat' % (count, plist[0], plist[1])
        np.savetxt(outname, pfreqwnp.reshape(20,20), fmt='%4f')
        count+=1
        '''
        # multiply by correlation values
        cedensity += (pfreqwnp * corrdict[k][0]) # ce score 
        dyndensity += (pfreqwnp * corrdict[k][1]) # df score
        freqdensity += pfreqwnp_o

    cedensityn = cp.normminmax(cedensity/freqdensity).reshape(20,20)
    dyndensityn = cp.normminmax(dyndensity/freqdensity).reshape(20,20)
    freqdensityn = cp.normminmax(freqdensity).reshape(20,20)

    np.savetxt(outprefix+'.freqmat', freqdensityn, fmt='%.2f')
    np.savetxt(outprefix+'.cemat', (cedensityn+cedensityn.T))
    np.savetxt(outprefix+'.dynmat', (dyndensityn+dyndensityn.T))
    cp._info('save %s.{cemat, dynmat, freqmat} with full precision' % outprefix)
    

# sum up all the .{cemat,dynmat} from aapaircorrdensity
# aa index: cp.aat01
def densitysum(args):
    assert len(args) == 2, 'Usage: pyton proc_minip.py densitysum filelist{.cemat/dynmat} outfile' 
    filelist = args[0]
    outfile= args[1]

    density = np.zeros((400,2))
    count=0
    for vfile in cp.loadlines(filelist):
        density += np.loadtxt(vfile)
        count+=1
        if(count%100==0):
            cp._info('%d families is done.' % count)

    np.savetxt(outfile, 1.0*(density)/np.max(density), fmt='%.4f')
    cp._info('save sum density to %s' % outfile)
    

# input single column data with cp.aas01 order
# put into square form with cp.aat01 order
# hard coded just for minip procedure
def densitysingle(args):
    assert len(args) == 2, 'Usage: python proc_minip.py densitysingle flat.vec out.mat'
    infile = args[0]
    outfile = args[1]

    flatnp = np.loadtxt(infile)
    aas01dict = dict(('%s %s' % (cp.aas01[i], cp.aas01[j]), flatnp[i*20+j]) for i in range(len(cp.aas01)) for j in range(len(cp.aas01)))
    aat01mat = np.array([aas01dict['%s %s' % (cp.aat01[i], cp.aat01[j])] for i in range(len(cp.aat01)) for j in range(len(cp.aat01))]).reshape(20,20)
    np.savetxt(outfile, 1.0*(aat01mat)/np.max(aat01mat), fmt='%.4f')
    cp._info('save aat01 mat to %s' % outfile)







def foo(args):
    print(len(args))

if __name__=='__main__':
    cp.dispatch(__name__)