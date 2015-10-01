# -*- coding: utf-8 -*-
"""
Created on Thu Oct 01 00:40:27 2015

@author: kjia
"""

# -*- coding: utf-8 -*-
from cgroup import cgroup
import itertools
import sys

def main():
    if len(sys.argv)< 3:
        print "Usage: proc_N-mer.py cluster_file n_merFile n_mer"
        return 
    infile = sys.argv[1]
    outfile = sys.argv[2]
    n_mer = int(sys.argv[3])
    print n_mer

    fo = open(outfile, 'w')
    with open(infile) as fp:
    	for line in fp:
    		cg = cgroup(line.strip())
    		print cg.getString()
    		if cg.getSize() == n_mer:
    			fo.write(cg.getString())
	    		print cg.getString()
	    	elif cg.getSize() < n_mer:
	    		continue
    		else:	
    			# generate combinations
    			for idx in list(itertools.combinations(range(cg.getSize()),n_mer)):
    				sub_cg = cgroup()
    				sub_cg.pdb = cg.pdb
    				for i in idx: # iterate all the tuples
    					sub_cg.AAgroup = sub_cg.AAgroup + cg.AAgroup[i]
    					sub_cg.resi.append(cg.resi[i])
    				fo.write(sub_cg.getString())
  	fp.close()
    fo.close()


if __name__=="__main__":
	main()    