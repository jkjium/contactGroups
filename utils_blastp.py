import commp as cp

# generate batch blastp script file
# iterate all query sequences, sseq.stub
# iterate all gap penalties, in commp gap gapdict
# iterate all available substitution matrices, smlist.sub
# read output string from file, outstring.sub
# arguments: dbname, evalue cutoff
# output file. seq.fa_matrix_gapopening_gapextend_dbname.out
def blastfarm(arglist):
    if len(arglist) < 6:
        cp._err('Usage: python utils_blasp.py blastfarm seq.stub sm.stub outstring.stub dbname evalue outfile')
        
    seqlistfile = arglist[0]
    seqlist = cp.loadlines(seqlistfile)

    smlistfile = arglist[1]
    smlist = cp.loadlines(smlistfile)

    outstringfile = arglist[2]
    outstrings = cp.loadlines(outstringfile)

    dbname = arglist[3]
    evalue = arglist[4]

    outfile = arglist[5]
    fout = open(outfile, 'w')

    blastroot = '/home/kjia/apps/ncbi-blast-2.2.31+-src/c++/built/'

    # matrix deterimines binpath, -matrix parameter, and gap penalties
    # blastp built-in matrices, gap penalties are different
    for sm in smlist:
        if sm in cp.gapdict:
            binpath = '/home/kjia/apps/ncbi-blast-2.2.31+-src/c++/built/b62.ReleaseMT/bin/blastp'
            smdata = sm
            gaplist = cp.gapdict[sm]
        else:
            binpath = '/home/kjia/apps/ncbi-blast-2.2.31+-src/c++/built/%s.ReleaseMT/bin/blastp' % sm
            smdata = 'BLOSUM62' 
            gaplist = cp.gapdict['BLOSUM62']
        
        for fa in seqlist:
            for g in gaplist:
                fout.write('%s -query %s -db $%s -outfmt "%s" -evalue %s -matrix %s -gapopen %d -gapextend %d -out %s_%s_%d_%d.out\n' % (binpath, fa, dbname, outstrings[0], evalue, smdata, g[0], g[1], fa, sm, g[0], g[1]))
    fout.close()
    cp._info('save batch script to %s' % outfile)

# created on 20200703 for LZTFL1_FYCO1_IFNAR2
# update version of blastfarm
# generate batch blastp script file
# use given binpath
# use given output format
# iterate all query sequences, sseq.stub
# iterate all gap penalties, in commp gap gapdict
# iterate all available substitution matrices, smlist.sub
# read output string from file, outstring.sub
# arguments: dbname, evalue cutoff
# output file. seq.fa_matrix_gapopening_gapextend_dbname.out
def blastfarmfmt(arglist):
    if len(arglist) < 6:
        cp._err('Usage: python utils_blasp.py blastfarm seq.stub sm.stub outstring.stub dbname evalue outfile')
        
    seqlistfile = arglist[0]
    seqlist = cp.loadlines(seqlistfile)

    smlistfile = arglist[1]
    smlist = cp.loadlines(smlistfile)

    outstringfile = arglist[2]
    outstrings = cp.loadlines(outstringfile)

    dbname = arglist[3]
    evalue = arglist[4]

    outfile = arglist[5]
    fout = open(outfile, 'w')

    blastroot = '/home/kjia/apps/ncbi-blast-2.2.31+-src/c++/built/'

    # matrix deterimines binpath, -matrix parameter, and gap penalties
    # blastp built-in matrices, gap penalties are different
    for sm in smlist:
        if sm in cp.gapdict:
            binpath = '/home/kjia/apps/ncbi-blast-2.2.31+-src/c++/built/b62.ReleaseMT/bin/blastp'
            smdata = sm
            gaplist = cp.gapdict[sm]
        else:
            binpath = '/home/kjia/apps/ncbi-blast-2.2.31+-src/c++/built/%s.ReleaseMT/bin/blastp' % sm
            smdata = 'BLOSUM62' 
            gaplist = cp.gapdict['BLOSUM62']
        
        for fa in seqlist:
            for g in gaplist:
                fout.write('%s -query %s -db $%s -outfmt "%s" -evalue %s -matrix %s -gapopen %d -gapextend %d -out %s_%s_%d_%d.out\n' % (binpath, fa, dbname, outstrings[0], evalue, smdata, g[0], g[1], fa, sm, g[0], g[1]))
    fout.close()
    cp._info('save batch script to %s' % outfile)

# main routine
if __name__ == '__main__':
	cp.dispatch(__name__)
