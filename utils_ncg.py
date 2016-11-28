import sys
from atom import atom
from protein import protein
from ncg import ncg

def writencg():
	if len(sys.argv) < 4:
		print 'writencg(): write non parametric contact group matrix for a (coarse-grained) pdb with size cutoff'
		print 'writencg(): python utils_ncg.py writencg 1t3r.pdb 3'
		return	

	pdbfile = sys.argv[2]
	size = int(sys.argv[3])
	outfile = pdbfile+'.ncg'

	print 'writencg(): pdbfile: %s' % pdbfile
	print 'writencg(): ncg size: %d' % size
	print 'writencg(): output: %s' % outfile

	ncgArray = []
	p = protein(pdbfile)
	for a in p.atoms:
		c = ncg(a, size)

def main():
	if len(sys.argv)<2:
		print 'Usage: python utils_ncg.py cmd pdbfile [args ...]'
		return

	dispatch = {
		'writencg': writencg
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
