import commp3 as cp
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