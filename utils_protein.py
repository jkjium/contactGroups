'''
get resid list of varname
'''
import sys
from protein import protein
from AAmap import AAmap

def resn2bfactor():
	if len(sys.argv) < 3:
		print 'resn2bfactor(): replace b factor values with residue type.'
		print 'resn2bfactor(): used for pymol spectrum b'
		return
	scoreValue = {
							'X':0,'-': 0,'.': 0,'A': 1,'C': 2,'D': 3,'E': 4,'F': 5,'G': 6,'H': 7,'I': 8,'K': 9,
							'L': 10,'M': 11,'N': 12,'P': 13,'Q': 14,'R': 15,'S': 16,'T': 17,'V': 18,'W': 19,'Y': 20, 'B': 3
						}
	aamap = AAmap()

	pdbfile = sys.argv[2]
	p = protein(pdbfile)
	outfile = '%s_rb.pdb' % pdbfile[:-4]
	fout = open(outfile, 'w')
	for a in p.atoms:
		newBFactor = scoreValue[aamap.getAAmap(a.resName)]
		print 'new b-factor: [%s : %s] -> %d' % (a.resName, aamap.getAAmap(a.resName), newBFactor)
		a.tempFactor = newBFactor
		fout.write(a.writeAtom())
	fout.close()
	print 'Output file: %s' % outfile


def main():
	if len(sys.argv)<3:
		print 'Usage: python utils_protein.py cmd pdbfile [args ...]'
		return

	dispatch = {
		'resn2bfactor': resn2bfactor
	}

	cmd = sys.argv[1]

	for key in dispatch:
		if key == cmd:
			dispatch[key]()


if __name__ == '__main__':
	main()
