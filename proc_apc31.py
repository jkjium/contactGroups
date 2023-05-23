import commp as cp
import numpy as np

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
    return cedict, idstub,colslis

def _apc3_21(cedict, idstub, collist):
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

def apc3_21(args):
    assert len(args) == 2, 'Usage: python proc_apc32.py apc32  '
    couplingfile = args[0]
    pass


def asc32(args):
    pass

if __name__=='__main__':
    cp.dispatch(__name__)