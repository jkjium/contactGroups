import commp as cp
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
	

if __name__ == '__main__':
    cp.dispatch(__name__)


