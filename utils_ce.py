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


# append shadow zscores to cefile
def zscore_outlierfilter(args):
    assert len(args) == 3, 'Usage: python utils_ce.py zscore_outlierfilter cefile cecolumn_index outfile'
    cefile = args[0]
    col = int(args[1])
    outfile = args[2]

    celines = cp.loadlines(cefile)
    celist = [float(list(line.split())[col]) for line in celines]

    cutoff, adj_m, adj_s = _outlierfilter(celist, 0) # cutoff is not used here
    cp._info('adj_m: %.8f, adj_s: %.8f' % (adj_m, adj_s))
    outlist = ['%s %.8f' % (celines[i], (celist[i]-adj_m)/adj_s) for i in xrange(0, len(celist))]

    with open(outfile, 'w') as fout:
        fout.write('%s\n' % '\n'.join(outlist))
    cp._info('save to %s' % outfile)



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
    assert len(args) == 6, 'Usage: python utils_ce.py keymap valuefile valuefile_key_columns {0,2,3} stubfile stubfile_key_columns {0,1,4} patch_value{no_patch} outfile'

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

# append all columns from keyvaluefile to targetfile according to key mapping
# no patch options, merge whatever matches
def keymerge_old(args):
    assert len(args) == 5, 'Usage: python utils_ce.py keymerge targetfile targetfile_key_columns {0,2,3} keyvaluefile keyvaluefile_key_columns {0,1,4} outfile'

    targetfile = args[0]
    vkclist = [int(i) for i in args[1].split(',')]
    stubfile = args[2]
    skclist = [int(i) for i in args[3].split(',')]
    outfile = args[4]

    # load value dictionary
    value_dict = {}
    values = cp.loadlines(stubfile)
    for line in values:
        sarr = line.split(' ')
        key = ' '.join([sarr[i] for i in vkclist])
        value = ' '.join([sarr[j] for j in xrange(0,len(sarr)) if j not in vkclist]) # save values with key columns
        value_dict[key] = value

    # append key-value to target line 
    # same order with the target file
    outlist = []
    for line in cp.loadlines(targetfile):
        sarr = line.split(' ')
        key = ' '.join([sarr[i] for i in skclist])
        if key in value_dict:
            outlist.append('%s %s' % (line, value_dict[key]))
        else:
            print('key: %s not found in %s' %(key, stubfile))

    with open(outfile, 'w') as fout:
        fout.write('%s\n' % '\n'.join(outlist))
    cp._info('save %d records to %s' % (len(outlist), outfile))

# paste all records from two file together according to the matched key
# updated version of the previous keymerge
# no patch options, merge whatever matches
def keymerge(args):
    assert len(args) == 5, 'Usage: python utils_ce.py keymerge tsvfile1 key_columns {0,2} tsvfile2 key_columns {1,3} outfile'
    infile1 = args[0]
    kclist1 = [int(i) for i in args[1].split(',')]
    infile2 = args[2]
    kclist2 = [int(i) for i in args[3].split(',')]
    outfile = args[4]

    valuedict1 ={}
    for line in cp.loadlines(infile1):
        sa = line.split(' ')
        k = ' '.join([sa[i] for i in kclist1])
        valuedict1[k] = line

    outdict ={}
    for line in cp.loadlines(infile2):
        sa = line.split(' ')
        k = ' '.join([sa[i] for i in kclist2])
        if k in valuedict1:
            outdict[k] = '%s %s' % (valuedict1[k], line)

    with open(outfile, 'w') as fout:
        fout.write('%s\n' % '\n'.join([outdict[k] for k in outdict]))
    cp._info('save %d merged records to %s' % (len(outdict), outfile))

'''
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

'''
RCW adjustment for pairwise correlations
rcw(i,j) = (sum_v(i)+sum_v(j) - MI(i,j)) / (n-1)
MI_rcw(i,j) = MI(i,j)/rcw(i,j)
# calcuate rcw with given stub order
'''
def _rcw(cedict, idpairstub, colslist, c):
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

    # column-wise average
    MIax = {}
    for i in range(0, len(colslist)):
        marginMI = 0.0
        for j in range(0, len(colslist)):
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
                # output spectrum
                #MIpdict[k] = (apc, MIp, p, q, MIax[colslist[p]], MIax[colslist[q]], avgMI)
                #print '%s %.8f %.8f %.8f\n' % (k, cedict[k], apc, MIp)
    return [MIpdict[k] for k in idpairstub]

# same as _apc0 but output full information
def _apc_full(cedict, idpairstub, colslist):
    avgMI = sum(cedict.itervalues()) / len(cedict)

    # column-wise average
    MIax = {}
    for i in range(0, len(colslist)):
        marginMI = 0.0
        for j in range(0, len(colslist)):
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
                #MIpdict[k] = (apc, MIp)
                # output spectrum
                MIpdict[k] = (apc, MIp, p, q, MIax[colslist[p]], MIax[colslist[q]], avgMI)
                #print '%s %.8f %.8f %.8f\n' % (k, cedict[k], apc, MIp)
    return [MIpdict[k] for k in idpairstub]


'''
calculate apc for triplet measures (ii, tc, sdii)
output apc3 spectrum {avgA avgB avgC total_avg}
'''
def _apc3c(cedict, idstub, colslist, c=0):
    if c==0:
        total_avg = sum(cedict.itervalues()) / len(cedict)
    else:
        total_avg = cp.adjmean(list(cedict.itervalues()), c)
    #nc = cp.ncr(len(colslist)-1, 2)
    # column-wise average
    vpos = collections.defaultdict(float)
    vlistdict = collections.defaultdict(list)
    for k in idstub:
        s = k.split(' ')
        v = cedict[k]
        '''
        vpos[s[0]]+=v
        vpos[s[1]]+=v
        vpos[s[2]]+=v
        '''
        vlistdict[s[0]].append(v)
        vlistdict[s[1]].append(v)
        vlistdict[s[2]].append(v)

    outlist = []
    for k in idstub:
        s = k.split(' ')
        #outlist.append((vpos[s[0]]/nc, vpos[s[1]]/nc, vpos[s[2]]/nc, total_avg))
        m1 = cp.adjmean(vlistdict[s[0]], c)
        m2 = cp.adjmean(vlistdict[s[1]], c)
        m3 = cp.adjmean(vlistdict[s[2]], c)
        outlist.append((m1,m2,m3,total_avg))
    return outlist

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
        

'''
Input: a space separted file, column id(s) for calculating RCW
Output: ssv file in the format of {rcw rcw_ce | rcw2 rcw_ce2 | ...}
'''
adj_func_dict = {'rcw': _rcw, 'apc': _apc, 'apc3': _apc3, 'apc_full': _apc_full}

def adjustment(args):
    assert len(args) == 5, 'Usage: python utils_ce.py adjustment infile.ssv stub_columns{0,1} apc ce_columns{2,3,5} outfile'

    # inputs
    infile = args[0]
    keyclist = [int(c) for c in args[1].split(',')] # which two columns are the indices
    #assert len(keyclist) == 2, 'illegal number of index columns for paired keys'
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
        outlist = adj_func(cedict, idpairstub, colslist) # returns a list of pair tuple (apc value, final MIp)
        outlists.append(outlist)

    #print '\n'.join([' '.join([' '.join(['%.8f' % v for v in outlist[i]]) for outlist in outlists]) for i in xrange(0, len(idpairstub))])
    with open(outfile, 'w') as fout:
        flatstr = '\n'.join([' '.join([' '.join(['%.8f' % v for v in outlist[i]]) for outlist in outlists]) for i in xrange(0, len(idpairstub))])
        fout.write('%s\n' % flatstr)
    cp._info('save {apc mip i j avgMI_i avgMI_j total_avg, ...} to %s' % outfile)

# copy adjustment 
# just for apc3
def apc3(args):
    assert len(args) == 5, 'Usage: python utils_ce.py apc3 infile.ssv stub_columns{0,1} ce_columns{2,3,5} adjmean_cutoff outfile'

    # inputs
    infile = args[0]
    keyclist = [int(c) for c in args[1].split(',')] # which two columns are the indices
    ceclist = [int(c) for c in args[2].split(',')] # target columns for the adjustment
    adjmean_cutoff = float(args[3])
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
        # won't improve accuracy and too slow ...
        #outlist = _apc3c(cedict, idpairstub, colslist, adjmean_cutoff) # returns a list of pair tuple (apc value, final MIp)
        outlist = _apc3(cedict, idpairstub, colslist) # returns a list of pair tuple (apc value, final MIp)
        outlists.append(outlist)

    #print '\n'.join([' '.join([' '.join(['%.8f' % v for v in outlist[i]]) for outlist in outlists]) for i in xrange(0, len(idpairstub))])
    with open(outfile, 'w') as fout:
        flatstr = '\n'.join([' '.join([' '.join(['%.8f' % v for v in outlist[i]]) for outlist in outlists]) for i in xrange(0, len(idpairstub))])
        fout.write('%s\n' % flatstr)
    #cp._info('save {rcw rcw_ce rcw2 rcw_ce2 ...} to %s' % outfile)
    cp._info('save {m1,m2,m3,total_m}, {...} to %s' % outfile)


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
    assert len(args) == 4, 'Usage: python utils_ce.py cflat2ccmat .cflatfile index.cols{0,1} cevalue.col{13} outprefix {.tick, .n.ccmat, .ccmat}'
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
            k1 = '%d %d' % (idlist[i], idlist[j]) 
            k2 = '%d %d' % (idlist[j], idlist[i]) 
            if k1 in cedict:
                v = cedict[k1]
            elif k2 in cedict:
                v = cedict[k2]
            else:
                v = 0.0
            ccmat[i,j] = ccmat[j,i] = v

    # output raw matrix file
    outccmatfile = '%s.ccmat' % outprefix
    np.savetxt(outccmatfile, ccmat, fmt='%.6f')
    # output normalized matrix file
    outnccmatfile = '%s.n.ccmat' % outprefix
    np.savetxt(outnccmatfile, cp.normminmax(ccmat), fmt='%.6f')
    # output tick file
    outtickfile = '%s.ccmat.tick' % outprefix
    np.savetxt(outtickfile, np.array(idlist), fmt='%i')

    cp._info('save to %s{.ccmat, .n.ccmat .ccmat.tick}' % outprefix)

def topsetsdii3(args):
    assert len(args) == 4, 'Usage: python utils_ce.py topsetsdii3 sdii3.file msailist.file threshold outfile'
    sdiifile = args[0]
    msaifile = args[1]
    th = float(args[2])
    outfile = args[3]

    # sdii3 file format c1-c3: index, c4: score
    # [kjia@lhb-ps4 ~/workspace/cadherin/cad/sdii3/stage] 2021-03-24 02:01:00
    # r1,r2,r3,sdii3_zscore
    outlist = []
    resilist = cp.loadlines(msaifile)
    for line in cp.loadlines(sdiifile):
        sarr = line.split(' ')
        r1 = sarr[0]
        r2 = sarr[1]
        r3 = sarr[2]
        sdii3 = float(sarr[3])
        if(r1 in resilist) and (r2 in resilist) and (r3 in resilist) and (sdii3>=th):
            outlist.append(line)

    with open(outfile, 'w') as fout:
        fout.write('%s\n' % '\n'.join(outlist))
    cp._info('save %d triplets in %s' % (len(outlist), outfile))

def topsetsdii2(args):
    assert len(args) == 4, 'Usage: python utils_ce.py topsetsdii2 dca.file msailist.file threshold outfile'
    sdiifile = args[0]
    msaifile = args[1]
    th = float(args[2])
    outfile = args[3]

    # sdii3 file format c1-c3: index, c4: score
    # [kjia@lhb-ps4 ~/workspace/cadherin/cad/sdii3/stage] 2021-03-24 02:01:00
    # r1,r2,r3,sdii3_zscore
    outlist = []
    resilist = cp.loadlines(msaifile)
    for line in cp.loadlines(sdiifile):
        sarr = line.split(' ')
        r1 = sarr[0]
        r2 = sarr[1]
        dca = float(sarr[2])
        if(r1 in resilist) and (r2 in resilist) and (dca>=th):
            outlist.append(line)

    with open(outfile, 'w') as fout:
        fout.write('%s\n' % '\n'.join(outlist))
    cp._info('save %d triplets in %s' % (len(outlist), outfile))






if __name__=='__main__':
    cp.dispatch(__name__)
