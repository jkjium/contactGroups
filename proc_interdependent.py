import commp as cp
import numpy as np
from scipy.stats import skew

# read file contains focused residues
# flag with len(id contained in focuslist), intra_: len(check_column)==len(id in focuslist)
def focusce(arglist):
    if len(arglist) < 4:
        cp._err('Usage: python pro_interdependent.py focusce cefile(.cflat) check_columns(0,1) focusvecfile(resi.vec) outfile')
    cefile = arglist[0]
    clist  = [int(i) for i in arglist[1].split(',')]
    # read from a .vec file
    focuslist = cp.loadlines(arglist[2])
    outfile = arglist[3]

    fout = open(outfile, 'w')
    for line in cp.loadlines(cefile):
        sarr = line.split(' ')
        count=0
        for i in clist:
            if sarr[i] in focuslist:
                count+=1
        fout.write('%d %s\n' % (count, line))
    fout.close()
    cp._info('save to %s' % outfile)




# general procedure for unskewing the ce distribution
# for pair, triplet, .....
# the last column is treated as the ce values
# rest of the columns compose of the key
# for example: 1 2 3 5.0 key: '1 2 3', value: 5.0
def unskew_zscore(arglist):
    if len(arglist) < 1:
        cp._err('Usage: python utils_interdependent.py unskew_zscore PF00003_full.txt.ce')
    cefile = arglist[0]

    keylist = []
    celist = []
    for line in cp.loadlines(cefile):
        sarr = line.split(' ')
        keylist.append('%s' % (' '.join(sarr[0:len(sarr)-1])))
        celist.append(float(sarr[len(sarr)-1]))
    oldskew = skew(celist)
    # take cubic root
    trans_ce = [x ** ( 1. / 3 ) for x in celist]
    newskew = skew(trans_ce)

    # calculate zscore
    nce = np.array(trans_ce)
    m = nce.mean()
    s = nce.std()
    zscore = [(x - m)/s for x in nce]

    outfile = '%s.trans.zscore' % (cefile)
    with open(outfile, 'w') as fout:
        for i in xrange(0, len(keylist)):
            fout.write('%s %.6f %.6f %.6f\n' % (keylist[i], celist[i], trans_ce[i], zscore[i]))
    print '%s %.6f %.6f %.6f %.6f' % (outfile, oldskew, newskew, m, s)


'''
for interdependent substitution project
based on the result of pfam31.0 2247 .cflat
'''
# input .dca file
#   $ head PF00003_full.txt.dca
# calculate the sknewness before and after taking cubic root transformation
# output .dca.zscore file based on the transformed values
#   format: r1 rn1 r2 rn2 dca transformed_dca transformed_zscore 
# print: skewness_before skewness_after
def unskew_zscore_dca(arglist):
    if len(arglist) < 1:
        cp._err('Usage: python utils_interdependent.py unskew_zscore_dca PF00003_full.txt.dca')
    dcafile = arglist[0]
    # r1 rn1 r2 rn2 mi dca
    # 57 A 58 A 0.615822 0.513832
    keylist = []
    dcalist = []
    for line in cp.loadlines(dcafile):
        sarr = line.split(' ')
        keylist.append('%s %s %s %s' % (sarr[0], sarr[1], sarr[2], sarr[3]))
        dcalist.append(float(sarr[5]))
    oldskew = skew(dcalist)
    # take cubic root
    trans_dca = [x ** ( 1. / 3 ) for x in dcalist]
    newskew = skew(trans_dca)

    # calculate zscore
    ndca = np.array(trans_dca)
    m = ndca.mean()
    s = ndca.std()
    zscore = [(x - m)/s for x in ndca]

    outfile = '%s.trans.zscore' % (dcafile)
    with open(outfile, 'w') as fout:
        for i in xrange(0, len(keylist)):
            fout.write('%s %.6f %.6f %.6f\n' % (keylist[i], dcalist[i], trans_dca[i], zscore[i]))
    print '%s %.6f %.6f %.6f %.6f' % (outfile, oldskew, newskew, m, s)

# input .mip file
# $ head PF00003_p90.mip
#   374 375 0.76824159 0.31576902 0.45247258
# output: .trans.mip.zscore (include mi and mip)
# format: 
def unskew_zscore_mip(arglist):
    if len(arglist) < 1:
        cp._err('Usage: python proc_interdependent.py unskew_zscore_mip PF00003_p90.mip')
    infile = arglist[0]
    keylist = []
    milist= []
    miplist =[]

    # p1  p2  mi   apc   mip
    # 374 375 0.76824159 0.31576902 0.45247258
    for line in cp.loadlines(infile):
        sarr = line.split(' ')
        keylist.append('%s %s' % (sarr[0], sarr[1]))
        milist.append(float(sarr[2]))
        miplist.append(float(sarr[4]))
    # save original skewness
    oldmiskew = skew(milist)
    oldmipskew = skew(miplist)

    # mi take cubic root
    #trans_mi = [ x ** ( 1. / 3 ) for x in milist ]
    trans_mi = [ cp.root(x,3) for x in milist ]
    newmiskew = skew(trans_mi)
    # mip take cubic root
    #trans_mip = [ x ** ( 1. / 3 ) for x in miplist ]
    trans_mip = [ cp.root(x,3) for x in miplist ]
    newmipskew = skew(trans_mip)

    # calculate zscore mi
    nmi = np.array(trans_mi)
    m_mi = nmi.mean()
    s_mi = nmi.std()
    zscore_mi = [(x-m_mi)/s_mi for x in nmi]
    # calculate zscore mip
    nmip = np.array(trans_mip)
    m_mip = nmip.mean()
    s_mip = nmip.std()
    zscore_mip = [(x-m_mip)/s_mip for x in nmip]

    # output mi
    outfile = '%s.trans.zscore' % (infile)
    with open(outfile, 'w') as fout:
        for i in xrange(0, len(keylist)):
            fout.write('%s %.6f %.6f %.6f %.6f %.6f %.6f\n' % (keylist[i], milist[i], trans_mi[i], zscore_mi[i], miplist[i], trans_mip[i], zscore_mip[i]))
    print '%s %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f' % (outfile, oldmiskew, newmiskew, m_mi, s_mi, oldmipskew, newmipskew, m_mip, s_mip)


# combine with previous .rcflat files
# adjust the order of values
# Format:
# pdbid chainID residueID1 residueID2 rn1 rn2 side_chain_center_distance ca_distance tip_atom_distance PfamID alignment_index1 alignment_index2 mip_zscore dca_zscore area1 area2
# use 'r1 r2' in trans_zscore as the key to match with .rcflat
# 3cyv A 5 6 K N 8.46205608 10.86056453 3.81297325 PF01208 146 147 -191 108.51798010
def unskew_dca_merge(arglist):
    if len(arglist) <3:
        cp._err('Usage: python proc_interdependent.py unskew_dca_merge original.rcflat .dca.trans.zscore output.newname.cflat')
    oldcflatfile = arglist[0]
    zscorefile = arglist[1]
    outfile = arglist[2]

    zscoredict = {}
    # format = 'r1 rn1 r2 rn2 dca trans_dca zscore'
    for line in cp.loadlines(zscorefile):
        sarr = line.split(' ')
        key = '%s %s' % (sarr[0], sarr[2])
        value = sarr[4:]
        zscoredict[key] = value
    
    with open(outfile, 'w') as fout:
        # pdbid chainID ri1 ri2 rn1 rn2 sc_dist ca_dist tip_dis PfamID p1 p2 mip_zscore dca_zscore area1 area2
        for line in cp.loadlines(oldcflatfile):
            sarr = line.split(' ')
            key = '%s %s' % (sarr[2], sarr[3])

            # check possible exceptions
            if key in zscoredict and sarr[12] == '-191':
                cp._err('key available but -191 found for: \n%s\n' % line)
            if (sarr[12] != '-191') and (key not in zscoredict):
                cp._err('Key not avaible but sarr[12] has value for: \n%s\n' % line)

            trans_dca_values = ''
            # keep -191 consistent
            if key in zscoredict:
                trans_dca_values = zscoredict[key]
            else:
                trans_dca_values == '-191 -191'
            #                                                                      0.pdbid  1.chain  2.ri1    3.ri2    4.rn1    5.rn2    6.sc_d   7.ca_d   8.tip_d  9.area1   10.area2 | 12.pfam  13.p1     14.p2     15.dca_z0 16.mip_z0 | 17.dca 18.dca_rcubic 19.dca_z1
            cflatstr = '%s %s %s %s %s %s %s %s %s %s %s | %s %s %s %s %s | %s' % (sarr[0], sarr[1], sarr[2], sarr[3], sarr[4], sarr[5], sarr[6], sarr[7], sarr[8], sarr[14], sarr[15],  sarr[9], sarr[10], sarr[11], sarr[13], sarr[12],   trans_dca_values)
            fout.write('%s\n' % (cflatstr))
    cp._info('save to %s' % outfile)


# read MSA.fa
# input target header (human P53)
# output the number of alphabet(amino acid) difference with the input sequence
def hamming_diff_with(arglist):
    if len(arglist) <2:
        cp._err('Usage: python proc_pf00870.pf hamming_diff_with msafile.fa target_header')
    msafile = arglist[0]
    target = arglist[1]
    targetseq = ''

    msalist = []
    for head, seq in cp.fasta_iter(msafile):
        if target in head:
            targetseq = seq
        msalist.append((head, seq))
    if targetseq == '':
        cp._err('target %s not found' % target)
    for head, seq in msalist:
        diff = [(str(i),targetseq[i],seq[i]) for i in xrange(0, len(targetseq)) if targetseq[i]!=seq[i]]
        diff_idx = ','.join(d[0] for d in diff)
        target_seg = ''.join([d[1] for d in diff])
        seq_seg = ''.join(d[2] for d in diff)
        print '%s %d %s %s %s' % (head, len(diff), target_seg, seq_seg, diff_idx)

if __name__=='__main__':
    cp.dispatch(__name__)