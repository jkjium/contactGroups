# -*- coding: utf-8 -*-
from protein import protein
# pre-process all atom pdb
# extract CA in Chain A and save to file
def main():
    fin = open('pdblist.txt', 'r')
    lines = fin.readlines()
    fin.close()
    
    for i in xrange(0,len(lines)):  
        line = lines[i].strip()
        pdb_filename=line+'.pdb'
        print pdb_filename
        p=protein(pdb_filename,'CA_A')
        p.writeChainACA('ca_'+line+'.pdb')

    pass
if __name__=="__main__":
	main()