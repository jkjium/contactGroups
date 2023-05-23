import commp as cp
from utils_pfammsa import pfammsa
import collections

def cm(args):
    assert len(args) == 3, 'Usage: python proc_cm.py cm PF00186.fa PF00186.vec19 PF00186.cm'
    fafile = args[0]
    cefile = args[1] # vec19
    outfile = args[2]

    cp._info('getting contacts ...')
    ctlist= []
    for t in cp.loadtuples(cefile):
        if float(t[18]) <= 4.5:
            ctlist.append((t[0],t[1],t[2],t[3],t[4],t[5])) # resi, msai, scorei
    cp._info('%d contacts found.' % len(ctlist))

    cp._info('extracting contact seqs ...')
    pfm = pfammsa(fafile)
    ctdict = collections.defaultdict(set)
    for i in range(len(pfm.msalist)):
        s = pfm.msalist[i][1]
        for c in ctlist:
            k = '%s,%s,%s,%s,%s,%s' % (c[0], c[1], c[2], c[3], c[4], c[5])
            a = s[int(c[2])]
            b = s[int(c[3])]
            if a!='.' and b!='.':
                ctdict[k].add((a,b)) 
    
    cp._info('collecting amino acid flips ...')
    outlist = []
    for k in ctdict:
        ctseqs = list(ctdict[k])
        for p in range(len(ctseqs)):
            for q in range(p+1, len(ctseqs)):
                a1 = ctseqs[p][0]
                a2 = ctseqs[p][1]
                b1 = ctseqs[q][0]
                b2 = ctseqs[q][1]
                if a1 == b2 and b1 == a2:
                    outlist.append('%s,%s%s,%s%s' % (k, a1,a2,b1,b2))
    
    with open(outfile, 'w') as fout:
        fout.write('%s\n' % '\n'.join(outlist))
    cp._info('save to %s' % outfile)

def cmseq(args):
    assert len(args) == 3, 'Usage: python proc_cm.py cmseq PF00186.aa.scoremat.fa 47,73,334,491,46,72,FA,AF PF00186.cm.FA.seq'
    scorefafile = args[0]
    cmstr = args[1]
    outfile = args[2]

    t = cmstr.split(',')
    i1 = int(t[4])
    i2 = int(t[5])

    n=1
    outlist = []
    for e in cp.fasta_iter(scorefafile):
        s = e[1]
        ctseq = '%s%s' % (s[i1], s[i2])
        if ctseq == t[6] or ctseq == t[7]:
            outlist.append('%d %s %s %s %s %s' % (n, s[i1-10:i1], s[i1], s[i1+1:i2], s[i2], s[i2+1:i2+10]))
        n+=1
    
    with open(outfile, 'w') as fout:
        fout.write('%s\n' % '\n'.join(outlist))
    cp._info('save to %s' % outfile)

if __name__=='__main__':
    cp.dispatch(__name__)