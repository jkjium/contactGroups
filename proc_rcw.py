import sys
import commp as cp
import numpy as np

'''
==> PF00098_p90.mip <== 
id1 id2 mi apc mip
7 8 0.02396475 0.00938778 0.01457697
7 9 0.10794275 0.08943567 0.01850708
7 10 0.07998329 0.09719031 -0.01720702
'''
# internal function
# updated from proc_mip2.py::mi2dict(sdii)
def _mi2dict(cefile):
    cedict={}
    idpairstub=[]
    idxset=set()
    for line in cp.loadlines(cefile):
        sarr = line.split(' ')

        k = '%s %s' % (sarr[0], sarr[1]) # use the first two column; sarr[0] < sarr[1] order
        idpairstub.append(k) # save the same order as the input file
        # save ce score (MI)
        cedict[k] = float(sarr[2]) 
        # get all column index in the infile
        idxset.add(int(sarr[0]))
        idxset.add(int(sarr[1]))

    colslist = list(idxset)
    colslist.sort()
    return cedict, idpairstub,colslist

def rcw(args):
    assert len(args) == 2
    infile = args[0]
    outfile = args[1]

    # load ce index scores from infile
    cedict, idpairstub, colslist = _mi2dict(infile)
    n=len(cedict)

    # calculate ith column MI sum
    sumdict = {}
    for i in xrange(0, len(colslist)):
        ith_column_sum = 0.0
        for j in xrange(0, len(colslist)):
            if i==j:
                continue
            k = ('%s %s' % (colslist[i], colslist[j])) if i<j else ('%s %s' % (colslist[j], colslist[i]))
            if k in cedict:
                ith_column_sum+=cedict[k]
        sumdict[colslist[i]]=ith_column_sum

    # calculate mi_rcw and write to outfile
    fout = open(outfile, 'w') 
    for k in idpairstub:
        col = [int(s) for s in k.split(' ')]
        # calculate rcw(i,j) = (v(i)+v(j) - M(i,j)) / (n-1)
        rcw = (sumdict[col[0]] + sumdict[col[1]] - cedict[k]) / (n-1)
        # MI_rcw(i,j) = MI(i,j)/rcw(i,j)
        mi_rcw = cedict[k] / rcw
        fout.write('%.8f %.8f\n' % (rcw, mi_rcw))

    fout.close()
    cp._info('save to %s' % outfile)

if __name__=='__main__':
    rcw(sys.argv[1:])