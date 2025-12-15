import commp as cp
import numpy as np
from collections import defaultdict

# input: the output from get_ortholog_by_id(args)
# 3-column tsv file, the last column may contain ','
# output flatten ',' column to {id, ortholog_id}
def flatten_orthologs(args):
    assert len(args) == 2, 'Usage: python proc_foxd.py flatten_orthologs 10.foxe.ad_ortho.tsv outfile'
    infile = args[0]
    outfile = args[1]

    outlist = []
    for t in cp.loadtuples(infile, delimiter='\t'):
        if t[2]=='_NA_':
            continue
        for s in t[2].split(','):
            outlist.append('%s\t%s\t%s' % (t[0],t[1],s))

    with open(outfile, 'w') as fout:
        fout.write('%s\n' % '\n'.join(outlist))
    cp._info('save flatten file to %s\n' % outfile)


# $ awk -F '\t' '{print NF}' nt__v__sp.tsv|sort -u = 3
def get_ortholog_by_id(args):
    assert len(args) == 3, 'Usage: python get_ortholog_by_id nt__v__ad.tsv 00.in.foxe.list.tsv outfile'
    orthofile = args[0]
    namefile =args[1]
    outfile = args[2]

    outdict = defaultdict(list)
    # nv_std_id alison_id
    for t in cp.loadtuples(namefile, delimiter='\t'):
        stdid = t[0]
        alisonid = t[1]
        for r in cp.loadtuples(orthofile, delimiter='\t'): # has 3 columns
            if alisonid in [s.strip() for s in r[1].split(',')]:
                outdict['%s\t%s' % (stdid, alisonid)].extend([s.strip() for s in r[2].split(',')])
        if len(outdict['%s\t%s' % (stdid, alisonid)]) == 0:
            outdict['%s\t%s' % (stdid, alisonid)].append('_NA_')
        
    with open(outfile, 'w') as fout:
        fout.write('%s\n' % '\n'.join(['%s\t%s' % (k, ','.join(outdict[k])) for k in sorted(outdict.keys())]))
    cp._info('save %d orthologs to %s' % (len(outdict), outfile))


# ortho one-to-one
# extract one-to-one orthologs from orthofinder outputs
# kjia@DESKTOP-K0I3VAM ~/workspace/foxd/orthofinder/Orthologues_ad 2024-12-15 18:31:00
# $ head ad__v__nt.tsv
# Orthogroup      ad      nt
# OG0000000       adig-s0038.g19  SVEP1-like-40
# OG0000000       adig-s0006.g329, adig-s0092.g75, adig-s0037.g103
# !! second column (ad) can have duplications, need to merge into a 1toMany list
# output: dict[ad] = nt
def _ortho_1to1(tuple_list):
    # merge potential 1toMany terms
    dict_1tomany = defaultdict(list)
    for t in tuple_list:
        query = t[1].split(',')
        hit = t[2].split(',')
        if(len(query)==1 and len(hit)==1):
            dict_1tomany[query[0]].append(hit[0])
    # remove 1 to many 
    return dict ((q, dict_1tomany[q][0]) for q in dict_1tomany if len(dict_1tomany[q])==1)

# load all orthofinder files and output 
def ortho_1to1_all(args):
    assert len(args) == 2, 'Usgae: python proc_foxd.py ortho_1to1_all stubfile outprefix'
    stubfile = args[0]
    outpref = args[1]
    ds = {}
    ds['names'] = cp.loadlines(stubfile)
    for n in ds['names']:
        ds[n] = _ortho_1to1(cp.loadtuples(n, '\t')[1:])
    # find common keys
    common_set = set(ds[ds['names'][0]].keys())
    for n in ds['names'][1:]:
        common_set &= set(ds[n].keys())
    cp._info('%d common keys found' % len(common_set))
    for n in ds['names']:
        outfile = '%s.%s.tsv' % (outpref, n) 
        with open(outfile , 'w') as fout:
            fout.write('%s\n' % '\n'.join(['%s\t%s' % (k, ds[n][k]) for k in common_set]))
        cp._info('save filtered 1to1 ortholog to %s' % outfile)

# ut _ortho_1to1()
def ortho_1to1(args):
    assert len(args) == 2, 'Usgae: python proc_foxd.py ortho_1to1 ad__v__nt.tsv outfile'
    infile = args[0]
    outfile = args[1]
    ortho_pair_list = cp.loadtuples(infile, '\t')[1:]
    out_1to1_dict = _ortho_1to1(ortho_pair_list)
    with open(outfile, 'w') as fout:
        fout.write('%s\n' % '\n'.join(['%s\t%s' % (k, out_1to1_dict[k]) for k in out_1to1_dict]))
    cp._info('save %d 1to1 ortho pair to %s' % (len(out_1to1_dict), outfile))


def foo(args):
    print(args)


if __name__=='__main__':
    cp.dispatch(__name__)
