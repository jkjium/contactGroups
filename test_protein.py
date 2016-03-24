'''
test protein class 
'''
from protein import protein


def main():
	pdbname = '1t3r.pdb'
	p = protein(pdbname)
	#p.writeCA(p.pdb+'.ca')
	p.printPDB()
	p.writeChainATips('AAtips.def', p.pdb+'.tip')

if __name__ == '__main__':
	main()
