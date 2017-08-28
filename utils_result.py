'''
post-process for result
'''
import sys
from protein import protein

# parse sdii_resi result
# ('B552', 'T')-('B604', 'S')-('B618', 'R') 1.198224093447873
#
def getresset():
	if len(sys.argv) < 2:
		print 'getresset(): python utils_result.py getresset result_sdii'
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
		'getresset': getresset
	}

	cmd = sys.argv[1]

	for key in dispatch:
		if key == cmd:
			dispatch[key]()


if __name__ == '__main__':
	main()
