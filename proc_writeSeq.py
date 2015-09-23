from protein import protein
# extract AA sequence and write to file as fasta format
def main():
    fin = open('pdblist.txt', 'r')
    lines = fin.readlines()
    fin.close()
    
    for i in xrange(0,len(lines)):  
        line = lines[i].strip()
#        pdb_filename=line+'.pdb'
        pdb_filename=line+'.tip'
        p=protein(pdb_filename,'fasta')
        p.writeSeq(line+'.fa')

    pass
if __name__=="__main__":
	main()# -*- coding: utf-8 -*-

