import itertools
import commp as cp

# for protsub testcase cath.h1
# filter out all the sequences included in the training pfam dataset
# input: cathS20.pfam.list
# output: cathS20.pfamid.list
def filterpfam(arglist):
    if len(arglist) < 3:
        cp._err('Usage: python proc_protsub.py filterpfam cathS20.pfam.list training.pfam.list cathS20.pfamid.list')
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


# generate .pool from cathS20.pfamfilter.list
# output: pdb pairs with cath family ID
# 1od2A01 1pixA02 3-90-226-10
def poolgen(arglist):
    if len(arglist) < 2:
        cp._err('Usage: python proc_protsub.py poolgen caths20.pfamfilter.list 3 caths20.pfamfilter.pair.list')
    infile = arglist[0]
    num = int(arglist[1]) # sequence number per family
    outfile =arglist[2]
    # [[12asA00 cathS20.00000.3-30-930-10.fa PF03590 3-30-930-10],[]]
    plist = [line.split(' ') for line in cp.loadlines(infile)]
    '''
    4-10-80-40
    =======================================
    ['2wdqA04', 'cathS20.07070.4-10-80-40.fa', 'na', '4-10-80-40']
    ['3gqhA02', 'cathS20.08941.4-10-80-40.fa', 'PF11962', '4-10-80-40']
    '''
    outlist = []
    for c, clist in itertools.groupby(plist, lambda x: (x[3])):
        '''
        print c
        print '======================================='
        print c, len(list(clist))
        '''
        count = num
        seqlist = list(clist)
        for i in xrange(0, len(seqlist)):
            for j in xrange(i+1, len(seqlist)):
                if count >0:
                    outlist.append('%s %s %s' % (seqlist[i][0], seqlist[j][0], seqlist[i][3]))
                    count-=1
    with open(outfile, 'w') as fout:
        fout.write('%s\n' % ('\n'.join(outlist)))
    cp._info('save to %s' % outfile)

if __name__ == '__main__':
	cp.dispatch(__name__)
