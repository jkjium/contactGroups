import commp as cp

# generate combination of pairs from a single column file
def pairgen(args):
    assert len(args) == 2, 'Usage: python proc_tuplegen.py pairgen infile.vec outfile.vec2'
    infile = args[0]
    outfile = args[1]
    ws = cp.loadlines(infile)
    outlist = ['%s %s' % (ws[i], ws[j]) for i in range(0, len(ws)) for j in range(i+1, len(ws))]
    with open(outfile, 'w') as fout:
        fout.write('%s\n' % '\n'.join(outlist))
    cp._info('save pairfile in %s' % outfile)

if __name__=='__main__':
    cp.dispatch(__name__)