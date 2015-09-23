# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 06:50:38 2015

@author: kjia
"""

# -*- coding: utf-8 -*-
from protein import protein
import sys

def main():
    if len(sys.argv)< 2:
        print "Usage: proc_SinglePDBFilter.py pdb(tip)file >> tip_clusters.txt"
        return 
  
    p=protein(sys.argv[1], 'v3',center='TIP')
    p.filterClusters()
if __name__=="__main__":
	main()    