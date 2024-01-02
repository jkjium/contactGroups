import commp3 as cp
import numpy as np

# bipartite apc procedure
# apc_rc = ((row_sum/col_n) outer (col_sum/row_n))/total_sum
# sm: numpy array score matrix; (pandas.to_numpy())
# r row by c column score matrix
# plt.figure(figsize=(8,6));plt.scatter(sm.flatten(), apc.flatten());plt.show()
def _apc_rc(sm):
    r,c = sm.shape
    '''
    row_mean = sm.sum(axis=1) / c # individual row mean
    col_mean = sm.sum(axis=0) / r # individual column mean
    apc = np.outer(row_mean, col_mean) / (sm.sum() / (r*c))
    '''
    return np.outer(sm.sum(axis=1)/c, sm.sum(axis=0)/r) / (sm.sum()/(r*c))


# call samap main procedure
# sn1,sn2: species names
# sf1,sf2: sam object file names
# map_p: homolog graph path
# resolutions: leiden clustering resolution for each species
# assignments: = {"bf":"cluster", "mm":"tissue_celltype"}
# sm = _run_samap_with_sam_obs('ad', '00.adsp.adig.counts_pr.h5ad', 'sp', '00.adsp.spis.counts_pr.h5ad', map_p='maps')
def _run_samap_with_sam_obs(sn1, sf1, sn2, sf2, map_p='maps/', resolutions=None, assignments=None):
    from samap.mapping import SAMAP
    from samap.analysis import (get_mapping_scores, GenePairFinder, sankey_plot, chord_plot, CellTypeTriangles, ParalogSubstitutions, FunctionalEnrichment, convert_eggnog_to_homologs, GeneTriangles)
    from samap.utils import save_samap, load_samap
    from samalg import SAM
    import pandas as pd

    cp._info('loading SAM objects.')
    sam1 = SAM()
    sam1.load_data(sf1)
    sam2 = SAM()
    sam2.load_data(sf2)
    sams = {sn1: sam1, sn2: sam2}

    cp._info('running SAMap ...')
    sm = SAMAP(sams, f_maps=map_p, resolutions=resolutions, keys=assignments)

    # neigh_from_keys = {'bf':True, 'mm':True}
    neigh_from_keys = dict((k,True) for k in assignments.keys()) if assignments!=None else None
    sm.run(pairwise=True, neigh_from_keys=neigh_from_keys)
    #save_samap(sm, 'ad_nv.samap.pkl')
    cp._info('done.')
    return sm




# !! incomplete
# convert prost output to blast fmt6 format
# adig-s0001.g1.t2        Hvul_1_1.g28296                         5523.5  1.12e-10
# ['adig-s0001.g1.t2', 'Hvul_1_1.g28296', '', '', '', '5523.5', '1.12e-10'] 
# blast -outfmt 6 # https://www.metagenomics.wiki/tools/blast/blastn-output-format-6
# 1      2      3      4      5        6       7      8     9     10    11    12
# qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
def prost2fmt6(args):
    assert len(args) == 2, 'Usage: python proc_coral_samap.py prost2fmt6 prost.out.tsv ad_to_hy.txt'
    infile = args[0]
    outfile = args[1]

    max_d = 0.0
    for t in cp.loadtuples(infile, delimiter='\t'):
        d = float(t[5])
        if d > max_d: max_d = d
    print(max_d)

    pr = np.genfromtxt(infile, dtype='str')
    print(pr)
    print(pr[:,2].astype(float))
    max_d = np.max(pr[:,2].astype(float))
    print(np.max(pr[:,2].astype(float)))

    pr[:,2]=np.power(-np.log(pr[:,3].astype(float)), 2)
    print(pr)
    # !! incomplete

##-------------------------------------------------------------------------------------
# gene_align_*: remove abnormal AA alphabets gaps = ['.','-',' ','*']
#             : reformat fasta headers; unique ID + transcription ID
##-------------------------------------------------------------------------------------
#
# kjia@DESKTOP-L0MMU09 ~/workspace/samap_coral/stage.nv_ad 2023-12-16 17:23:01
# $ awk '{print substr($1,1,5)}' nvec.genes.tsv|sort|uniq -c
#   33973 Nvec_
# change fasta header to align with nvec.genes.txt 
# from:">v1g152167" to: ">Nvec_v1g152167"
def gene_align_nvec(args):
    assert len(args) == 2, 'Usage: python proc_coral_samap.py gene_align_nvec xenSp1.proteins.fa xesp.proteome.rename.fa'
    infile = args[0]
    outfile = args[1]

    outlist = []
    for h,s in cp.fasta_iter(infile):
        ts = s.translate(str.maketrans('','',''.join(cp.gaps))) # remove '.' at the end of some sequences
        # h = '>v1g152167'
        th = 'Nvec_%s' % h
        outlist.append('>%s\n%s' % (th,ts))
    
    with open(outfile, 'w') as fout:
        fout.write('%s\n' % '\n'.join(outlist))
    cp._info('save to %s' % outfile)


# kjia@DESKTOP-L0MMU09 ~/workspace/samap_coral/stage.xe_ad 2023-12-15 15:19:35
# $ awk '{print substr($1,1,4)}' *genes.tsv|sort|uniq -c
#  29015 Xesp
#   1795 orph
# Xesp records have the same length
# $ grep Xesp *genes.tsv|awk '{print length($1)}'|sort |uniq -c
#  29015 12
# kjia@DESKTOP-L0MMU09 ~/workspace/samap_coral/stage.xe_ad 2023-12-15 15:21:08
# $ grep ">" xenSp1.proteins.fa | awk '{print length($1)}'|sort|uniq -c
#  28010 23

# change fasta header to align with adig.rownames.txt 
# from:">Xe_029129-T2 Xe_029129" to: "Xesp_029129.t2"
def gene_align_xesp(args):
    assert len(args) == 2, 'Usage: python proc_coral_samap.py gene_align_xesp xenSp1.proteins.fa xesp.proteome.rename.fa'
    infile = args[0]
    outfile = args[1]

    import re
    outlist = []
    for h,s in cp.fasta_iter(infile):
        ts = s.translate(str.maketrans('','',''.join(cp.gaps))) # remove '.' at the end of some sequences
        # h = '>Xe_029168-T1 Xe_029168'
        sarr = re.split(r'[_\-\s]',h) # ['Xe', '029168', 'T1', 'Xe', '029168']
        #print(sarr)
        th = 'Xesp_%s.%s' % (sarr[1], sarr[2].lower())
        outlist.append('>%s\n%s' % (th,ts))
    
    with open(outfile, 'w') as fout:
        fout.write('%s\n' % '\n'.join(outlist))
    cp._info('save to %s' % outfile)


# change fasta header to align with adig.rownames.txt 
# from:">Sc4wPfr_165.g10029.t1" to: "Hvul_g10029_1.t1"
def gene_align_hvul(args):
    assert len(args) == 2, 'Usage: python proc_coral_samap.py gene_align_hvul hydra2.0_genemodels.aa hvul.proteome.rename.fasta'
    infile = args[0]
    outfile = args[1]

    outlist = []
    for h,s in cp.fasta_iter(infile):
        ts = s.translate(str.maketrans('','',''.join(cp.gaps))) # remove '.' at the end of some sequences
        #print(">%s\n%s\n" % (h,ts))
        sarr0 = h.split('.') # ['>Sc4wPfr_165', 'g10029', 't1']
        th = 'Hvul_%s_1.%s' % (sarr0[1], sarr0[2]) # Hvul_g10029_1.t1
        outlist.append('>%s\n%s' % (th,ts))
    
    with open(outfile, 'w') as fout:
        fout.write('%s\n' % '\n'.join(outlist))
    cp._info('save to %s' % outfile)


# change fasta header to align with adig.rownames.txt 
# from:">adig_s0001.g1.t2 gene=adig_s0001.g1" to: "adig-s0001.g1.t2"
def gene_align_adig(args):
    assert len(args) == 2, 'Usage: python proc_coral_samap.py gene_align_adig adig.pep.fasta output.fastas'
    infile = args[0]
    outfile = args[1]

    outlist = []
    for h,s in cp.fasta_iter(infile):
        ts = s.translate(str.maketrans('','',''.join(cp.gaps))) # remove '.' at the end of some sequences
        #print(">%s\n%s\n" % (h,ts))
        sarr0 = h.split(' ') # ['>adig_s0001.g1.t2' 'gene=adig_s0001.g1']
        sarr1 = sarr0[0].split('_') # ['adig', 's0001.g1.t2']
        th = '%s-%s' % (sarr1[0], sarr1[1]) # adig-s0001.g1.t2
        outlist.append('>%s\n%s' % (th,ts))
    
    with open(outfile, 'w') as fout:
        fout.write('%s\n' % '\n'.join(outlist))
    cp._info('save to %s' % outfile)



def foo(args):
    print(args)

if __name__=='__main__':
    cp.dispatch(__name__)