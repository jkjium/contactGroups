import commp as cp
import numpy as np
from protein import protein
from atom import atom

def filterpdbwithresi(args):
	assert len(args)==3, 'Usage: python proc_ompt.py filterpdbwithresi pdbfile resi.vec outfile'
	pdbfile = args[0]
	resilist = [int(str_resi) for str_resi in cp.loadlines(args[1])]
	outfile = args[2]

	outlist =[]
	p = protein(pdbfile)
	for at in p.atoms:
		if at.resSeq in resilist:
			outlist.append(at.writeAtom())

	with open(outfile, 'w') as fout:
		fout.write('\n'.join(outlist))

	cp._info('save %d atoms in filtered pdb to %s' % (len(outlist), outfile))

def combineadjmat(args):
	assert len(args)==3, 'Usage: python proc_ompt.py combineadjmat {dca.adjmat,contact.adjmat} {0.5,1.0} outfile'
	adjfilelist = args[0].split(',')
	valuelist = [float(v) for v in args[1].split(',')]
	outfile = args[2]
	if len(adjfilelist)!= len(valuelist):
		cp._err('length mismatch: adjmat list: %d, value list: %d ' % (len(adjfilelist), len(valuelist)))
	adjmatlist = [np.loadtxt(fn) for fn in adjfilelist]
	for i in range(len(valuelist)):
		adjmat = adjmatlist[i]
		adjmat[adjmat==1.0] = valuelist[i]

	totalsum = np.sum(adjmatlist, 0)	
	np.savetxt(outfile, totalsum, fmt='%.4f')
	cp._info('save sum to %s' % outfile)

if __name__ == '__main__':
    cp.dispatch(__name__)


