import commp as cp
import numpy as np
import collections

# calculate cutoff for informative elements
# introduced in sdii shadow
# return cutoff, adjusted mean and std (for z-score calculation)
def _outlierfilter(nlist, threshold):
    nplist = np.array(nlist)
    outlier = nplist.mean()+nplist.std()

    nplist_no_outlier = np.array([v for v in nlist if v < outlier])
    m = nplist_no_outlier.mean()
    s = nplist_no_outlier.std()
    #cutoff = nplist_no_outlier.mean() + threshold * nplist_no_outlier.std()
    cutoff = m + threshold * s
    return cutoff, m, s


# filter top ce values using outlier procedure introduced in sdii shadow
# updated version of proc_topsdii_stuple.py
# output selected ce lines with adjusted z-score appended at the last row (instead of just column tuples)
# works for both pairs and triplets
def topce_outlierfilter(args):
    assert len(args)==4, '\n\nUsage: python utils_ce.py topce_outlierfilter cefile{PF00870.cflat} cecolumn_index{6} outlier_threshold{3} outfile'
    # inputs
    cefile = args[0]
    col = int(args[1]) # 0-based
    threshold = float(args[2])
    outfile = args[3]

    # A space separated value file
    celines = cp.loadlines(cefile)
    celist = [float(list(line.split(' '))[col]) for line in celines]

    cutoff, adj_m, adj_s = _outlierfilter(celist, threshold)
    # output celines with a ce value larger than the cutoff
    outlist = ['%s %.8f' % (celines[i], (celist[i]-adj_m)/adj_s) for i in xrange(0, len(celist)) if celist[i] > cutoff]

    with open(outfile, 'w') as fout:
        fout.write('%s\n' % '\n'.join(outlist))
    cp._info('%s: cutoff: %.4f #ofIPV: %d/%d' % (outfile, cutoff, len(outlist), len(celines)))


# caculate the count of unique alphabets for given columns
# input: 
#   a csv file; 
#   a set of column indices, separated by ','. such as {1,2,3}
#   .stub file: output indices
# options: patch or output indices with a given stub (will be the xticks in the plot)
# output:
#   csv file in the order of the given .stub file 
#   format: resi msai count
# works for both pairs and triplets
def alphabetfreq(args):
    assert len(args)==5, '\n\nUsage: python utils_ce.py alphabetfreq PF00870.cflat 0,1,2 PF00870.xtick.stub {or na} mapfile {or na} outfile'

    infile = args[0]
    clist = [int(i) for i in args[1].split(',')]
    stubfile = args[2]
    mapfile = args[3]
    outfile = args[4]

    # get all elements from given columns
    elist = []
    for line in cp.loadlines(infile):
        sarr = line.split(' ')
        for i in clist:
            elist.append(sarr[i])

    # get unique alphabet counts with given stub
    countdict = collections.Counter(elist)
    stublist = [alphabet for alphabet in cp.loadlines(stubfile)] if stubfile != 'na' else countdict.keys()

    resimap = collections.defaultdict(lambda: '-191') 
    if mapfile!='na':
        for line in cp.loadlines(mapfile):
            sarr = line.split(' ')
            resimap[sarr[2]] = sarr[0]
        outlist = ['%s %s %d' % (resimap[k], k, countdict[k]) if k in countdict else '%s %s 0' % (resimap[k], k) for k in stublist]
    else: # no resi information
        outlist = ['%s %d' % (k, countdict[k]) if k in countdict else '%s 0' % (k) for k in stublist]

    with open(outfile, 'w') as fout:
        fout.write('%s\n' % '\n'.join(outlist))
    cp._info('save counts to %s' % outfile)

# general key mapping function, value filtered by stub with mapping the corresponding columns in two files
# given a stub file and columns {0,2,3} that compose of a key {space separated key}
# given columns {0,3,4} in the value file that map the key
# get all the values of line[key]
# option: patch -191 if value is less than stub, default: 'no_patch'
def keymap(args):
    assert len(args) == 6, 'Usage: python utils_ce.py keymap valuefile valuefile_key_columns {0,2,3} stubfile stubfile_key_columns {0,1,4} patch_value outfile'

    valuefile = args[0]
    vkclist = [int(i) for i in args[1].split(',')]
    stubfile = args[2]
    skclist = [int(i) for i in args[3].split(',')]
    patch_value = args[4]
    outfile = args[5]

    # load value dictionary
    value_dict = {}
    values = cp.loadlines(valuefile)
    for line in values:
        sarr = line.split(' ')
        key = ' '.join([sarr[i] for i in vkclist])
        value = ' '.join([sarr[j] for j in xrange(0,len(sarr)) if j not in vkclist])
        value_dict[key] = value
    num_v = len(values[0].split(' '))

    # stub key list
    stub_key_list = []
    for line in cp.loadlines(stubfile):
        sarr = line.split(' ')
        key = ' '.join([sarr[i] for i in skclist])
        stub_key_list.append(key)

    if patch_value == 'no_patch': # output min(value_key_list, mapped(stub_key_list)) 
        outlist = ['%s %s' % (k, value_dict[k]) for k in stub_key_list if k in value_dict] # only output available value, igore non-existing key
    else:
        outlist = ['%s %s' % (k, value_dict[k]) if k in value_dict else '%s %s' % (k, ' '.join([patch_value]*(num_v-len(skclist)))  ) for k in stub_key_list]
    
    with open(outfile, 'w') as fout:
        fout.write('%s\n' % '\n'.join(outlist))

    cp._info('save %d/%d records in %s' % (len(outlist), len(stub_key_list), outfile))

def _ce2dict(celines, keyclist, ceclist):
    cedict={}
    idpairstub=[]
    idxset=set()
    for line in celines:
        sarr = line.split(' ')

        k = ' '.join(['%s' % sarr[c] for c in keyclist]) # keyclist should always have two columns {id1, id2}
        # save the same order as the input file
        # this will control the order of output
        idpairstub.append(k) 

        # save ce scores 
        #print k, repr(ceclist), len(sarr)
        cedict[k] = [float(sarr[c]) for c in ceclist]

        # get all column index in the infile
        for c in keyclist:
            idxset.add(int(sarr[int(c)]))

    colslist = list(idxset)
    colslist.sort()
    return cedict, idpairstub,colslist

'''
RCW adjustment for pairwise correlations
rcw(i,j) = (sum_v(i)+sum_v(j) - MI(i,j)) / (n-1)
MI_rcw(i,j) = MI(i,j)/rcw(i,j)
# calcuate rcw with given stub order
'''
def _rcw(cedict, idpairstub, colslist):
    # calculate column-wise sum
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
    outlist = []
    for k in idpairstub:
        col = [int(s) for s in k.split(' ')]
        # calculate rcw(i,j) = (v(i)+v(j) - M(i,j)) / (n-1)
        rcw = (sumdict[col[0]] + sumdict[col[1]] - cedict[k]) / (len(cedict)-1)
        # MI_rcw(i,j) = MI(i,j)/rcw(i,j)
        mi_rcw = cedict[k] / rcw
        outlist.append((rcw, mi_rcw))    
    return outlist


'''
calculate MIp based on .rcol file and raw sdii(MI)
raw sdii format:
masi1 msai2 sdii
APC(a,b) = (MIax * MIbx) / avgMI
output: single column MIp corresponding to the original order in SDII file (for paste to append)
'''
def _apc(cedict, idpairstub, colslist):
    avgMI = sum(cedict.itervalues()) / len(cedict)

    # column-wise sum
    MIax = {}
    for i in xrange(0, len(colslist)):
        marginMI = 0.0
        for j in xrange(0, len(colslist)):
            if i == j:
                continue
            k = ('%s %s' % (colslist[i], colslist[j])) if i < j else ('%s %s' % (colslist[j], colslist[i]))
            if k in cedict:
                marginMI+= cedict[k]
        MIax[colslist[i]] = marginMI / (len(colslist)-1)
    
    # calculate MIp
    MIpdict={}
    MIp = 0.0
    MIplist=[]
    for p in xrange(0, len(colslist)):
        for q in xrange(p+1, len(colslist)):
            apc = (MIax[colslist[p]] * MIax[colslist[q]])/avgMI
            k = '%d %d' % (colslist[p], colslist[q])
            if k in cedict:
                MIp = cedict[k] - apc
                MIpdict[k] = (apc, MIp)
                #print '%s %.8f %.8f %.8f\n' % (k, cedict[k], apc, MIp)
    return [MIpdict[k] for k in idpairstub]


'''
Input: a space separted file, column id(s) for calculating RCW
Output: ssv file in the format of {rcw rcw_ce | rcw2 rcw_ce2 | ...}
'''
adj_func_dict = {'rcw': _rcw, 'apc': _apc}

def adjustment(args):
    assert len(args) == 5, 'Usage: python utils_ce.py adjustment infile.ssv stub_columns{0,1} ce_columns{2,3,5} rcw outfile'

    # inputs
    infile = args[0]
    keyclist = [int(c) for c in args[1].split(',')] # which two columns are the indices
    assert len(keyclist) == 2, 'illegal number of index columns for paired keys'
    adj_func = adj_func_dict[args[2]]
    ceclist = [int(c) for c in args[3].split(',')] # target columns for the adjustment
    outfile = args[4]

    # load ce index scores from infile
    # cedicts['2 4'] = [1.1, 2.2, 3.3]
    # idpairstub: contains the output order
    # colslist: all the indices in an ordered list
    cedicts, idpairstub, colslist = _ce2dict(cp.loadlines(infile), keyclist, ceclist)

    # calculate the adjustment one by one
    outlists = []
    for i in xrange(0, len(ceclist)):
        cedict = dict((k, cedicts[k][i]) for k in cedicts.keys())
        outlist = adj_func(cedict, idpairstub, colslist) # returns a list of pair tuple (rcw_value, ce_rcw)
        outlists.append(outlist)

    #print '\n'.join([' '.join([' '.join(['%.8f' % v for v in outlist[i]]) for outlist in outlists]) for i in xrange(0, len(idpairstub))])
    with open(outfile, 'w') as fout:
        flatstr = '\n'.join([' '.join([' '.join(['%.8f' % v for v in outlist[i]]) for outlist in outlists]) for i in xrange(0, len(idpairstub))])
        fout.write('%s\n' % flatstr)
    cp._info('save {rcw rcw_ce rcw2 rcw_ce2 ...} to %s' % outfile)


# generate ccmatrix for downstream analysis {EC analysis}
# input:
# index_colids: columns that gives tuple id, {resi or msai}, works for pair and triplets ..
# value_colid: the one column that contains ce value
# outprefix: two output files are generated, ccmat in np.loadtxt() format and the IDs {resi or msai} to label the matrix
# example:
# 
# kjia@kjia-PC ~/workspace/src 2020-08-24 17:05:01
# $ python utils_ce.py cflat2ccmat t.cflat 2,3 6 t.out
# 2020-08-24 17:05:13|41238|0|INFO|save to t.out {.ccmat, .ccmat.tick}
# t.cflat:
# -191 -191 459 460 H A 0.400267 19.225526 0.736970 7.552786 0.459538 0.013046 0.771686 0.181111
# -191 -191 459 461 H E 0.165445 7.425152 0.548973 4.681673 0.263010 -0.889776 0.640704 -0.789504
# -191 -191 460 461 A E 0.336884 16.040377 0.695814 6.924250 0.320346 -0.626382 0.684237 -0.466913
# t.out.ccmat:
# 0.0000 0.4003 0.1654
# 0.4003 0.0000 0.3369
# 0.1654 0.3369 0.0000
# t.out.ccmat.tick:
# 459
# 460
# 461
def cflat2ccmat(args):
    assert len(args) == 4, 'Usage: python utils_ce2ccmat .cflatfile index.cols{0,1} cevalue.col{13} outprefix {.cecol, .ccmat}'
    cflatfile = args[0]
    index_colids = [int(i) for i in args[1].split(',')]
    value_colid = int(args[2])
    outprefix = args[3]

    lines = cp.loadlines(cflatfile)

    # get index set
    ids=set()
    for line in lines:
        sarr = line.split()
        ids.update([int(sarr[i]) for i in index_colids])
    idlist=list(ids)
    idlist.sort()

    # get key: ids - value: target ce dictionary
    def _func_getkvpair(sarr, index_collids, value_colid):
        k = ' '.join([sarr[i] for i in index_colids])
        v = float(sarr[value_colid])
        return k,v
    cedict = dict(_func_getkvpair(line.split(), index_colids, value_colid) for line in lines)

    # fill in ccmat
    ccmat = np.zeros((len(idlist), len(idlist)))
    for i in range(len(idlist)):
        for j in range(i+1, len(idlist)):
            ccmat[i,j] = ccmat[j,i] = cedict['%d %d' % (idlist[i], idlist[j])]

    # output files
    outccmatfile = '%s.ccmat' % outprefix
    np.savetxt(outccmatfile, ccmat, fmt='%.4f')
    outtickfile = '%s.ccmat.tick' % outprefix
    np.savetxt(outtickfile, np.array(idlist), fmt='%i')

    cp._info('save to %s{.ccmat, .ccmat.tick}' % outprefix)



if __name__=='__main__':
    cp.dispatch(__name__)