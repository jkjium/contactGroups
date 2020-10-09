import commp as cp
import numpy as np 

def foo(args):
    cp._info(args)

# append dca,mi,mip zscores from PF01694-dca-mi-mip.cflat2 to 3q2v_ab.msai.resi.dist
def appendce(args):
    assert len(args) == 4, 'Usage: python proc_cadherin.py appendce 3q2v_ab.msai.resi.dist PF00028.cflat2 cecolumns{0,2,4} outfile'

    distfile = args[0]
    cflat2file = args[1]
    cols = [int(i) for i in args[2].split(',')] # 16,19,22 - zscore dca, mi, mip
    outfile = args[3]

    #npcflat2 = np.loadtxt(cflat2file, dtype=np.object, delimiter=' ')
    #print(npcflat2[:,cols])

    # load 3q2v_ab.msai.resi.dist
    # 504 338 A 79 H B 58 G 22.1298
    #_msaikey = lambda x: '%s %s' % (x[0], x[1]) 
    #distdict = dict((_msaikey(line.split()), line) for line in cp.loadlines(distfile))
    distdict = {}
    for line in cp.loadlines(distfile):
        sarr=line.split()
        keystr = '%s %s' % (sarr[3], sarr[6])
        if keystr in distdict:
            print('%s: %s already exist' % (keystr, distdict[keystr]) )
            cp._err('exit')
        else:
            distdict[keystr] = line
    cp._info('%d dist lines loaded.' % len(distdict))
    
    # load cflat2 
    cflat2dict = {}
    for line in cp.loadlines(cflat2file):
        sarr = line.split()
        valuestr = ' '.join([sarr[i] for i in cols])
        keystr1 = '%s %s' % (sarr[2], sarr[3])
        keystr2 = '%s %s' % (sarr[3], sarr[2])
        cflat2dict[keystr1] = valuestr
        cflat2dict[keystr2] = valuestr
    cp._info('%d cevalues loaded.' % len(cflat2dict))

    # appenc values by key mapping
    outstrlist = []
    for k in distdict:
        appendvalue = cflat2dict[k] if k in cflat2dict else '-191 -191 -191'
        outstrlist.append('%s %s' % (distdict[k], appendvalue))
    cp._info('%d values appended' % len(outstrlist))

    # output new .dist file 
    with open(outfile, 'w') as fout:
        fout.write('%s\n' % '\n'.join(outstrlist))
    cp._info('save appended dist file to %s' % outfile)


def vec32mat(args):
    assert len(args) == 2, 'Usage: python proc_cadherin.py vec32mat input.vec3 outprefix'
    vec3file = args[0]
    outprefix = args[1]

    ticks=set()
    vdict={}
    for line in cp.loadlines(vec3file):
        sarr = line.split()
        ticks.add(int(sarr[0]))
        ticks.add(int(sarr[1]))
        keystr = '%d %d' % (int(sarr[0]), int(sarr[1]))
        vdict[keystr] = int(sarr[2]) # 0 or 1

    ticklist = list(ticks)
    ticklist.sort()
    print ticklist

    mat = []
    for i in range(len(ticklist)):
        r=[]
        for j in range(len(ticklist)):
            k = '%d %d' % (ticklist[i], ticklist[j])
            r.append(vdict[k])
        mat.append(r)
    npmat = np.array(mat)

    outticks = '%s.ticks' % outprefix
    with open(outticks, 'w') as fout:
        fout.write('%s\n' % ('\n'.join([str(i) for i in ticklist])))
    cp._info('save ticks to %s' % outticks)

    outmat = '%s.mat' % outprefix
    np.savetxt(outmat, npmat, fmt='%i')
    cp._info('save matrix to %s' % outmat)


# split sequences (without fasta header) into separate fa files
def split(args):
    fasfile = args[0]
    count = 0
    for line in cp.loadlines(fasfile):
        head = '3q2v65_%04d' % (count)
        outfile = head+'.fa'
        with open(outfile, 'w') as fout:
            fout.write('>%s\n%s\n' % (head, line))
        count+=1
    cp._info('%d fa files saved.' % count)


# input.anlflat: output from utils_testcase.py testpool
# resilist: 0-based index number of target residues
# /home/kjia/workspace/others/sayane/hydrophobic/stage/3q2v65.pool_needle_EBLOSUM62_10.0_0.5_.alnflat
def alnflatresi(args):
    assert len(args) == 2, 'Usage: python proc_cadherin.py alnflatresi input.anlflat{output from testpool} 0,1,22,24'
    alnfile = args[0]
    resilist = [int(c) for c in args[1].split(',')]
    # palign
    from utils_testcase import palign
    totalset = set()
    for line in cp.loadlines(alnfile):
        p = palign(line)
        idx = 0
        outstrlist = []
        for i in range(len(p.seqA)):
            if p.seqA[i] != '-':
                if idx in resilist:
                    outstrlist.append(p.seqB[i])
                    #print('%s.' % p.seqB[i])
                    #print('%d %s.' % (i, p.seqA[i]))
                idx+=1
        outstr=' '.join(outstrlist)
        #print("\n")
        print(outstr)

if __name__=='__main__':
    cp.dispatch(__name__)