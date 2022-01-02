import commp as cp
import numpy as np
from scipy.stats import rankdata
import string
import itertools
from protein import protein

# column selection & stat
# input: .vec19, .map (for stat) ec_cutoff(percentage rank)
# output: .scols file, .scols.stat
def scols(args):
    assert len(args) == 5, 'Usage: python utils_ps3.py scols .vec19 .map ec_col_id ec_cutoff outprefix'
    vecfile = args[0]
    mapfile = args[1]
    ecid = int(args[2]) # 0-based
    cutoff = float(args[3])
    outprefix = args[4]

    # load map
    m = cp.loadlines(mapfile)

    outlist = []
    msaiset = set()
    zscorelist = []
    for v in cp.loadtuples(vecfile):
        if float(v[ecid])<=cutoff:
            outlist.append('%s %s %s %s %s %s' % (v[2], v[3], v[4], v[5], v[ecid-1], v[ecid])) # m1 m2 i1 i2 zscore p-rank
            msaiset.add(v[2])
            msaiset.add(v[3])
            zscorelist.append(float(v[ecid-1])) # get all zscores
    
    # scol output
    outscolfile = '%s.scols' % outprefix
    with open(outscolfile, 'w') as fout:
        fout.write('%s\n' % '\n'.join(outlist))

    # stat info
    npzscore = np.array(zscorelist)
    statstr = '%s %.2f %d %d %.4f %.4f %.4f' % (outprefix[:7], cutoff, len(m), len(msaiset), 1.0*len(msaiset)/len(m), npzscore.mean(), npzscore.std())
    outstatfile = '%s.stat' % outscolfile
    with open(outstatfile, 'w') as fout:
        fout.write('%s\n' % statstr)
    cp._info('save to {%s, %s}, stat %s' % (outscolfile, outstatfile, statstr))


# append rand percentile to each ec score
# output a vec19 file
def appendrankp(args):
    assert len(args) == 2, 'Usage: python utils_ps3.py trans_rankp .vec15 out.vec15'
    infile = args[0]
    outfile = args[1]
    vec15 = np.loadtxt(infile)

    n = vec15.shape[0]
    # just hard coded ...
    mi = vec15[:,6] 
    dca = vec15[:,8] 
    mip = vec15[:,10] 
    dcap = vec15[:,12]

    miz = vec15[:,7]
    dcaz = vec15[:,9] 
    mipz = vec15[:,11] 
    dcapz = vec15[:,13]

    #mir = 1.0*((-mi).argsort()+1)#/n
    #mir = ((-mi).argsort()).argsort()+1#/n
    mir = rankdata(-mi)/n
    #print(np.stack((mi,-mi, mir),axis=1))
    dcar = rankdata(-dca)/n
    mipr = rankdata(-mip)/n
    dcapr = rankdata(-dcap)/n

    indices = vec15[:,0:6]
    dist = vec15[:,14]

    # put back together
    outce = np.stack((mi,miz,mir,dca,dcaz,dcar,mip,mipz,mipr,dcap,dcapz,dcapr,dist), axis=1)
    outnp = np.concatenate((indices, outce), axis=1)

    np.savetxt(outfile, outnp, fmt='%d %d %d %d %d %d %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f')
    cp._info('save to %s' % outfile)



# for protsub testcase cath.h1
# filter out all the sequences included in the training pfam dataset
# input: cathS{20}.pfam.list
# output: cathS{20}.pfamid_label.list, append label to previous list, 1: included in the training dataset; 2: not included
#   17gsA01.fa.json PF02798 1
def filterpfam(arglist):
    if len(arglist) < 3:
        cp._err('Usage: python utils_protsub.py filterpfam cathS20.pfam.list training.pfam.list cathS20.pfamid.list')
    infile = arglist[0]
    filterfile = arglist[1]
    outfile = arglist[2]

    filterset = set(cp.loadlines(filterfile))
    cp._info('%d PfamID in the training dataset' % len(filterset))

    outlist = []
    # cathS20.00000.3-30-930-10.fa.json PF03590
    for line in cp.loadlines(infile):
        sarr = line.split(' ')
        pfset = set(sarr[1].split(','))
        flag = len(pfset.intersection(filterset))
        outlist.append('%s %d' % (line, flag))

    with open(outfile ,'w') as fout:
        fout.write('%s\n' % ('\n'.join(outlist)))
    cp._info('save to %s' % outfile)


# generate .pool from caths{40,60,95}.filtered.pdb-cathid.list
# output: pdb pairs with cath family ID
#   1od2A01 1pixA02 3-90-226-10
def poolgen(arglist):
    if len(arglist) < 2:
        cp._err('Usage: python proc_protsub.py poolgen caths20.pfamfilter.list 3 caths20.pfamfilter.pair.list')
    infile = arglist[0]
    num = int(arglist[1]) # sequence number per family
    outfile =arglist[2]
    # [[12asA00 3-30-930-10],[...]]
    plist = [line.split(' ') for line in cp.loadlines(infile)]
    outlist = []
    for c, clist in itertools.groupby(plist, lambda x: (x[3])):
        count = num
        seqlist = list(clist)
        for i in xrange(0, len(seqlist)):
            for j in xrange(i+1, len(seqlist)):
                if count >0:
                    outlist.append('%s %s %s' % (seqlist[i][0], seqlist[j][0], seqlist[i][1]))
                    count-=1
    with open(outfile, 'w') as fout:
        fout.write('%s\n' % ('\n'.join(outlist)))
    cp._info('save to %s' % outfile)


# for extracting cath sequences {S40, S60, only}
# extract sequences by list of pdb IDs
# output: one file for each sequence in fasta format
#   for pfamscan to filter out sequences contain traning pfam domains
def extractseq(args):
    assert len(args) == 2, 'Usage: python utils_ps3.py extractseq pdbs_s40-20.list cath-dataset-nonredundant-S40.atom.fa'
    pdblistfile = args[0]
    fafile = args[1]

    # load pdblist
    pdblist = cp.loadlines(pdblistfile)
    # filter ambiguous alphabet
    #trans = string.maketrans(''.join(cp.ambaa), ''.join(['.' for i in xrange(len(cp.ambaa))]))
    count = 0
    for head, seq in cp.fasta_iter(fafile):
        sarr = head.split('/')
        pdbid = sarr[0][-7:]
        if pdbid in pdblist:
            with open(pdbid+'.fa', 'w') as fout:
                #fout.write('>%s\n%s\n' % (head, seq.translate(trans).upper()))
                fout.write('>%s\n%s\n' % (head, seq.upper()))
            count+=1
            if count % 1000 == 0:
                cp._info('%d seq saved.' % count)
    cp._info('[%d] .fa saved.' % count)


# the result of tmalign is two structures in one pdb file
# to use the old rmsd pipeline they need to be separated into single structure files
# input: aln.1c9rA01.1mu2B01.pdb 
# output: aln.1c9rA01.1mu2B01.A.pdb aln.1c9rA01.1mu2B01.B.pdb
def splitalnab(args):
    assert len(args) == 1, 'Usage: python utils_ps3.py splitalnab aln.1c9rA01.1mu2B01.pdb'
    alnpdbfile = args[0]
    pa = protein(alnpdbfile, chain="A")
    if len(pa.atoms)==0:
        cp._err('%s: no atoms in chain A' % alnpdbfile)
    pb = protein(alnpdbfile, chain="B")
    if len(pb.atoms)==0:
        cp._err('%s: no atoms in chain B' % alnpdbfile)
    
    pa.writepdb(alnpdbfile[:-4]+".A.pdb")
    cp._info('save chain A to %s' % alnpdbfile[:-4]+".A.pdb")
    pb.writepdb(alnpdbfile[:-4]+".B.pdb")
    cp._info('save chain B to %s' % alnpdbfile[:-4]+".B.pdb")

    

                
                
if __name__ == '__main__':
        cp.dispatch(__name__)
