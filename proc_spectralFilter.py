# -*- coding: utf-8 -*-
"""
Created on Thu Sep 10 03:43:38 2015

@author: kjia
"""

# -*- coding: utf-8 -*-
from protein import protein
import sys

def main():
    if len(sys.argv)< 3:
        print "Usage: proc_spectralFilter.py pdb(tip)file cluster_size>> tip_clusters.txt"
        return 
    p=protein(sys.argv[1], 'v4',center='TIP')
    n = int(sys.argv[2])
    p.spectralClustering(n)
if __name__=="__main__":
	main()    