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


def pdbcut():
	if len(sys.argv) < 5:
		print 'pdbcut(): write pdb by residue segment'
		print 'pdbcut(): python utils_protein.py pdbcut 1t3r.pdb A 5-15'
		print 'pdbcut(): python utils_protein.py pdbcut 1t3r.pdb all 5-15'
		return

	pdbfile = sys.argv[2]
	chain = sys.argv[3]
	rangeStr = sys.argv[4]

	rangeArray = rangeStr.split('-')

	rBegin = int(rangeArray[0])
	rEnd = int(rangeArray[1])

	pdbname = pdbfile[0:4]
	outfile = '%s_%s_%d_%d.rpdb' % (pdbname, chain, rBegin, rEnd)

	print 'pdbcut():pdbfile: %s' % pdbfile
	print 'pdbcut():pdb: %s' % pdbname
	print 'pdbcut():chain: %s' % chain
	print 'pdbcut():residue range: %d - %d' % (rBegin, rEnd)

	p = protein(pdbfile)
	out = []
	if chain == 'all':
		for a in p.atoms:
			if a.resSeq <= rEnd and a.resSeq >= rBegin:
				out.append(a)
	else:
		for a in p.atoms:
			if (a.resSeq <= rEnd and a.resSeq >= rBegin and a.chainID == chain):
				out.append(a)

	fout = open(outfile, 'w')
	print 'pdbcut():output: %s' % outfile
	print 'pdbcut():%d atoms written.' % len(out)
	for a in out:
		fout.write(a.writeAtom())
	fout.close()

def main():
	if len(sys.argv)<2:
		print 'Usage: python utils_protein.py cmd pdbfile [args ...]'
		return

	dispatch = {
		'resn2bfactor': resn2bfactor, 'pdbcut': pdbcut
	}

	cmd = sys.argv[1]

	if sys.argv[1] not in dispatch:
		print 'invalid cmd: %s' % sys.argv[1]
		return

	for key in dispatch:
		if key == cmd:
			dispatch[key]()


if __name__ == '__main__':
	main()
