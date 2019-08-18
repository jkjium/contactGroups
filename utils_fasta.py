import commp as cp
'''
input: namelist: list of fa names
output: save each .fa with its name
'''
def splitbyname(arglist):
    if len(arglist) < 2:
        cp._err('Usage: python utils_fasta.py splitbyname namelist full_fa_seq.txt')
    namefile = arglist[0]
    fastafile = arglist[1]

    # load name lines
    namelist = cp.loadlines(namefile)

    # load fasta
    fadict = dict()
    for h,s in cp.fasta_iter(fastafile):
        sarr = h.split(' ')
        fadict[sarr[0]] = s
        #print '%s\n%s' % (sarr[0],s)

    for name in namelist:
        outfile = '%s.fa' % (name)
        if name not in fadict:
            cp._info('target: %s not in dict' % name)
            continue
        with open(outfile, 'w') as fout:
            fout.write('>%s\n%s\n' % (name, fadict[name]))
    
if __name__ == '__main__':
        cp.dispatch(__name__)