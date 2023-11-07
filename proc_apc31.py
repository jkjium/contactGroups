import commp as cp
import numpy as np
import collections

# works for multivarable
# input: tsv lines, key column indices, ce value column indecies
# output: cedict['0 2 3']= (ii, tc, sdii), idstub = ['0 1 2', '0 1 3', ...,], collist = [0,1,2,3,...]
def _ce2dict(celines, keyclist, ceclist):
    cedict={}
    idstub=[]
    idxset=set()
    for line in celines:
        sarr = line.split(' ')

        k = ' '.join(['%s' % sarr[c] for c in keyclist]) # keyclist should always have two columns {id1, id2}
        # save the same order as the input file
        # this will control the order of output
        idstub.append(k) 

        # save ce scores 
        #print k, repr(ceclist), len(sarr)
        cedict[k] = [float(sarr[c]) for c in ceclist]

        # get all column index in the infile
        for c in keyclist:
            idxset.add(int(sarr[int(c)]))

    colslist = list(idxset)
    colslist.sort()
    return cedict, idstub,colslist

# the orginal apc3 without adjusted mean
def _apc3(cedict, idstub, colslist):
    total_avg = sum(cedict.itervalues()) / len(cedict)
    nc = cp.ncr(len(colslist)-1, 2)
    # column-wise average
    vpos = collections.defaultdict(float)
    for k in idstub:
        s = k.split(' ')
        v = cedict[k]
        vpos[s[0]]+=v
        vpos[s[1]]+=v
        vpos[s[2]]+=v
    outlist=[]
    for k in idstub:
        s = k.split(' ')
        outlist.append((vpos[s[0]]/nc, vpos[s[1]]/nc, vpos[s[2]]/nc, total_avg))
    return outlist

def apc3(args):
    assert len(args) == 4, 'Usage: python utils_ce.py apc3 PF00001.i3.ii.sdii.tc.dtc.mdist.vec8 0,1,2 5 tc.m4'

    # inputs
    infile = args[0]
    keyclist = [int(c) for c in args[1].split(',')] # which two columns are the indices
    ceclist = [int(c) for c in args[2].split(',')] # target columns for the adjustment
    outfile = args[3]

    # load ce index scores from infile
    # cedicts['2 4'] = [1.1, 2.2, 3.3]
    # idpairstub: contains the output order
    # colslist: all the indices in an ordered list
    cedicts, idstub, colslist = _ce2dict(cp.loadlines(infile), keyclist, ceclist)

    # calculate the adjustment one by one
    outlists = []
    cedict = dict((k, cedicts[k][0]) for k in cedicts.keys()) # get single cescore
    out4 = _apc3(cedict, idstub, colslist) # returns a list of pair tuple (m1,m2,m3,all_avg)
    outlists.append(out4)

    with open(outfile, 'w') as fout:
        flatstr = '\n'.join([' '.join([' '.join(['%.8f' % v for v in outlist[i]]) for outlist in outlists]) for i in xrange(0, len(idstub))])
        fout.write('%s\n' % flatstr)
    cp._info('save {m1,m2,m3,total_m} to %s' % outfile)


def _apc32(cedict, idstub, collist):
    total_avg = sum(cedict.itervalues()) / len(cedict)
    nc = len(collist)

    vpos = collections.defaultdict(float)
    for k in idstub: # all position triplets
        s = k.split(' ')
        v = cedict[k]
        vpos['%s %s' % (s[0], s[1])]+=v
        vpos['%s %s' % (s[0], s[2])]+=v
        vpos['%s %s' % (s[1], s[2])]+=v

    outlist = []
    for k in idstub:
        s= k.split(' ')
        k1 = '%s %s' % (s[0], s[1])
        k2 = '%s %s' % (s[0], s[2])
        k3 = '%s %s' % (s[1], s[2])
        outlist.append((vpos[k1]/nc, vpos[k2]/nc, vpos[k3]/nc, total_avg))
    return outlist


# only works on single cescore
def apc32(args):
    assert len(args) == 4, 'Usage: python utils_ce.py apc3 PF00001.i3.ii.sdii.tc.dtc.mdist.vec8 0,1,2 5 tc.m4'

    # inputs
    infile = args[0]
    keyclist = [int(c) for c in args[1].split(',')] # which two columns are the indices
    ceclist = [int(c) for c in args[2].split(',')] # target columns for the adjustment
    outfile = args[3]

    # load ce index scores from infile
    # cedicts['2 4'] = [1.1, 2.2, 3.3]
    # idpairstub: contains the output order
    # colslist: all the indices in an ordered list
    cedicts, idstub, colslist = _ce2dict(cp.loadlines(infile), keyclist, ceclist)

    # calculate the adjustment one by one
    outlists = []
    cedict = dict((k, cedicts[k][0]) for k in cedicts.keys()) # get single cescore
    out4 = _apc32(cedict, idstub, colslist) # returns a list of pair tuple (m1,m2,m3,all_avg)
    outlists.append(out4)

    with open(outfile, 'w') as fout:
        flatstr = '\n'.join([' '.join([' '.join(['%.8f' % v for v in outlist[i]]) for outlist in outlists]) for i in xrange(0, len(idstub))])
        fout.write('%s\n' % flatstr)
    cp._info('save {m1,m2,m3,total_m} to %s' % outfile)


if __name__=='__main__':
    cp.dispatch(__name__)