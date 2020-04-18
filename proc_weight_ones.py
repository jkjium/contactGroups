import commp as cp

# generate weight file with all values of one
# for pfam31.0 p90 stage .90.weight
def onesweight(arglist):
    if len(arglist) < 2:
        cp._info('Usage: python onesweight template_weight_file outfile')
    with open(arglist[1],'w') as fout:
        fout.write('%s\n' % '\n'.join(['1.00000000' for line in cp.loadlines(arglist[0])]))
    cp._info('save to %s' % arglist[1])
if __name__ == '__main__':
	cp.dispatch(__name__)