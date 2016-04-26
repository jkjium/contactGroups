'''
get resid list of varname
'''
import sys
from protein import protein

def getseq():
	if len(sys.argv) < 3:
		print 'getseq(): python utils_protein.py getseq pdbfile chainID'
		return

	pdbfile = sys.argv[2]
	chain = sys.argv[3]
	p = protein(pdbfile, chain)
	print p.seq

def main():
	if len(sys.argv)<3:
		print 'Usage: python utils_protein.py cmd pdbfile [args ...]'
		return

	dispatch = {
		'getseq': getseq
	}

	cmd = sys.argv[1]

	for key in dispatch:
		if key == cmd:
			dispatch[key]()


if __name__ == '__main__':
	main()
