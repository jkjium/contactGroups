# -*- coding: utf-8 -*-
import sys
from protein import protein
# pre-process all atom pdb
# extract CA in Chain A and save to file
def main():
    if len(sys.argv) < 3:
        print 'Usage: proc_getTip.py pdblist'
        return
    pdblist = sys.argv[1]
    fin = open(pdblist, 'r')
    lines = fin.readlines()
    fin.close()
    
    for i in xrange(0,len(lines)):  
        line = lines[i].strip()
        #pdb_filename=line+'.pdb'
        pdb_filename=line
        print pdb_filename
        p=protein(pdb_filename,'TIP')
        p.writeChainATips('AAtips.def',line+'.tip')

    pass
if __name__=="__main__":
	main()