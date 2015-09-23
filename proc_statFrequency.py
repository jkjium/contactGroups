#!/usr/bin/python
import MySQLdb as mdb
import sys

def main():
    clusterDict = {}    
    try:
#        con = mdb.connect('localhost', 'cluster', 'cluster', 'db_cluster_1', charset='utf8');
        con = mdb.connect('localhost', 'root', '', 'db_cluster_test', charset='utf8');
        cur = mdb.cursors.SSCursor(con)
#        cur = con.cursor()
        
        cur.execute('select alignstr, weight from t_clusters')

        count=1
        while True:
            row = cur.fetchone()
            if not row:
                break              
#            print "alignstr: [%s], weight: %f\n" % (row[0], row[1])
       
            if count%100000==0:
                print "%d clusters loaded\n" % count
            if row[0] in clusterDict:
                clusterDict[row[0]]+=row[1]
            else:
                clusterDict[row[0]]=row[1]
            count+=1
            
        


        fd=open('clusters.freq', 'w')
        for key in clusterDict:
            output= "%s %f\n" % (key, clusterDict[key]) 
            fd.write(output)
        fd.close()   
            
        print "write %d unique clusters\n" % len(clusterDict)

        
#        con.commit()
#        print "done."
#        cur.execute("SELECT VERSION()")
#    
#        ver = cur.fetchone()
        
#        print "Database version : %s " % ver
        
    except mdb.Error, e:
      
        print "Error %d: %s" % (e.args[0],e.args[1])
        sys.exit(1)
        
    finally:    

        if cur:
            cur.close()
            
        if con:    
            con.close()

    pass
if __name__=="__main__":
	main()
    