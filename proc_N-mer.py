# -*- coding: utf-8 -*-
"""
Created on Thu Oct 01 00:40:27 2015

@author: kjia
"""

# -*- coding: utf-8 -*-
from cgroup import cgroup
from protein import protein
import itertools
import sys

def main():
    if len(sys.argv)< 3:
        print "Usage: proc_N-mer.py tip_file cluster_file n_mer dist_cutoff >> n_mer_outfile"
        return 
    infile = sys.argv[1]
    infile2 = sys.argv[2]
    n_mer = int(sys.argv[3])
    cutoff = float(sys.argv[4])

    p = protein(infile, center='TIP')
    p.initCGResiMap()

    with open(infile2) as fp:
    	for line in fp:
    		cg = cgroup(line.strip())
    		if cg.getSize() == n_mer:
   				if p.cgResiGroupFilter(cg, cutoff) == True:
 						print cg.getString()
	    	elif cg.getSize() < n_mer:
	    		continue
    		else:	
    			# generate combinations
    			for idx in list(itertools.combinations(range(cg.getSize()),n_mer)):
    				sub_cg = cgroup()
    				sub_cg.pdb = cg.pdb
    				sub_cg.chain = cg.chain
    				for i in idx: # iterate all the tuples
    					sub_cg.AAgroup = sub_cg.AAgroup + cg.AAgroup[i]
    					sub_cg.resi.append(cg.resi[i])
    				if p.cgResiGroupFilter(sub_cg, cutoff) == True:
 						print sub_cg.getString()
  	fp.close()


if __name__=="__main__":
	main()    