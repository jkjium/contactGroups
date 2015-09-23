# -*- coding: utf-8 -*-
import MySQLdb as mdb
from protein import protein
from msa import msa
import sys
# pre-process all atom pdb
# extract CA in Chain A and save to file
def main():
    if len(sys.argv)< 4:
        print "Usage: proc_MySQL.py msa_file topology center"
        return

    msafile = sys.argv[1]
    
    topology = sys.argv[2]  
    c = sys.argv[3] 
    pdbfile = (sys.argv[1])[0:4]+'.'+c
    
    print '%s %s %s' % (msafile, topology, c)
    try:
        con = mdb.connect('localhost', 'cluster', 'cluster', 'db_cluster_1', unix_socket='/home/kjia/workspace/data/mysql/mysql.sock', charset='utf8');
        #con = mdb.connect('localhost', 'root', '', 'db_cluster_test', charset='utf8');
        cur = con.cursor()
        
#        p=protein('1alu.pdb', 'alpha')
        p=protein(pdbfile, topology, center=c)
        p.filterClusters()        
        
#        m=msa('1alu_PF00489.fa', p)
        m=msa(msafile, p)
        m.writeAlignClustersdb(cur)
        m.printName()
        con.commit()
#        print "done."
#        cur.execute("SELECT VERSION()")
#    
#        ver = cur.fetchone()
        
#        print "Database version : %s " % ver
        
    except mdb.Error, e:
      
        print "Error %d: %s" % (e.args[0],e.args[1])
        sys.exit(1)
        
    finally:    
            
        if con:    
            con.close()
if __name__=="__main__":
	main()
