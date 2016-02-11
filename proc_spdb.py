# generate concised pdb from tip(domain) file
import sys
from protein import protein
from atom import atom
def main():
	if len(sys.argv) < 2:
		print 'Usage python proc_spdb.py pdblist.txt' 
		return
	pdblist = sys.argv[1]
	fin = open(pdblist, 'r')
	lines = fin.readlines()
	fin.close()

	for i in xrange(0, len(lines)):
		line = lines[i].strip()
		tip_file = line + '.tip'
		print tip_file
		p = protein(tip_file, 'alpha', center='TIP')
		p.writeSPDB()

'''
	pdbfile = sys.argv[1]
	p = protein(pdbfile, 'alpha', center="TIP")
	p.writeSPDB()
'''

if __name__=="__main__":
	main()

