import commp as cp
import collections


# for sh2-btk interaction prediction
# combine MSAs 
def msapadding(args):
    if len(args)!= 4:
        cp._err('Usage: python proc_sh2btk.py msapadding from.fa to.fa tax.list outfile')

    def _gettax(header):
        # 3BP2_HUMAN/458-538
        sa = header.split('/')
        # 3BP2_HUMAN
        sa1 = sa[0].split('_')
        return sa1[1]

    fromfile = args[0]
    tofile = args[1]
    taxfile = args[2]
    outfile = args[3]

    seqbytaxlistdict = collections.defaultdict(list)
    seqbytaxnumdict = collections.defaultdict(int)

    # load common tax list
    commontaxlist = set(cp.loadlines(taxfile))

    # separate from.fa into different species
    for header, seq in cp.fasta_iter(fromfile):
        tax = _gettax(header)
        if tax in commontaxlist:
            seqbytaxlistdict[tax].append((header, seq))
    '''
    print len(seqbytaxlistdict)
    for k in seqbytaxlistdict:
        print k, len(seqbytaxlistdict[k])
    '''

    fout = open(outfile, 'w')
    # padding onto to.fa
    for header, seq in cp.fasta_iter(tofile):
        tax = _gettax(header)
        if tax in commontaxlist:
            index = seqbytaxnumdict[tax] % len(seqbytaxlistdict[tax])
            print tax, index
            paddingheader, paddingseq = seqbytaxlistdict[tax][index]
            outstr = '>%s|%s\n%s%s\n' % (header, paddingheader, seq, paddingseq)
            fout.write(outstr)
            seqbytaxnumdict[tax]+=1

    fout.close()
    cp._info('save to %s' %  outfile)


# for two msa shares the same header (gen + tax)
# pfam full alignments
def msapaddingbyheader(args):
    if len(args)!=4:
        cp._err('Usage: python proc_sh2btk.py msapaddingbyheader msa1.fa msa2.fa sharedgentax.list outfile')
    
    msa1 = args[0]
    msa2 = args[1]
    headerfile = args[2]
    outfile = args[3]

    seqbyheaderlist1 = collections.defaultdict(list)
    seqbyheaderlist2 = collections.defaultdict(list)

    commonheader = set(cp.loadlines(headerfile))
    # separate msa1
    for header, seq in cp.fasta_iter(msa1):
        # 3BP2_HUMAN/458-538
        sarr = header.split('/')
        subheader = sarr[0]
        if subheader in commonheader:
            seqbyheaderlist1[subheader].append((header, seq))
    
    # separate msa2
    for header, seq in cp.fasta_iter(msa2):
        sarr = header.split('/')
        subheader = sarr[0]
        if subheader in commonheader:
            seqbyheaderlist2[subheader].append((header, seq))


    fout = open(outfile, 'w')
    count = 0
    # generate combined MSA with all-to-all scheme for each common header
    for subheader in commonheader:
        for h1,s1 in seqbyheaderlist1[subheader]:
            for h2,s2 in seqbyheaderlist2[subheader]:
                fout.write('>%s|%s\n%s%s\n' % (h1,h2,s1,s2))
                count+=1
    fout.close()
    cp._info('save %d combined sequences in %s' % (count, outfile))




















if __name__=='__main__':
    cp.dispatch(__name__)