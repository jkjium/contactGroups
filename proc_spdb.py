# generate concised pdb from tip(domain) file
import sys
from protein import protein
from atom import atom
def main():
	if len(sys.argv) < 2:
		print 'Usage python proc_spdb.py XXXX_A.domain'
		return
	pdbfile = sys.argv[1]
	p = protein(pdbfile, 'alpha', center="TIP")
	p.writeSPDB()

if __name__=="__main__":
	main()

