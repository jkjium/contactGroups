import sys
import os
import time
import math
from itertools import groupby
import operator as op

import numpy as np
'''
from scipy.cluster.hierarchy import dendrogram
from scipy.cluster.hierarchy import linkage
from scipy.cluster.hierarchy import fcluster
'''

import itertools
import inspect
import collections
# for testing drun dcall
import subprocess
from functools import reduce

illab =['/', '+', '?', '*', '&', '$', '\\', ' ', '|']

gaps = ['.','-',' ','*']

# abnormal AA
abaa = ['.','-','*','?','X','Z','U','B','O','J']

# ambiguous AA
ambaa = ['X','Z','U','B','O','J','-']

aafull = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X', '-',
                 'a', 'r', 'n', 'd', 'c', 'q', 'e', 'g', 'h', 'i', 'l', 'k', 'm', 'f', 'p', 's', 't', 'w', 'y', 'v', 'b', 'z', 'x', '.']

# AA alphabet for emboss format
aa201 = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
                'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

# AA alphabet sorted by singlet name
aas01 = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
# append gap '.' into the list
aas02 = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '.']

# AA alphabet sorted by type
aat01 = ['A','I','L','V','M','F','W','G','P','C','N','Q','S','T','Y','D','E','R','H','K']

aa2a={'ARG':'R','HIS':'H','LYS':'K','ASP':'D','GLU':'E',
      'SER':'S','THR':'T','ASN':'N','GLN':'Q','CYS':'C',
          'SEC':'U','GLY':'G','PRO':'P','ALA':'A','VAL':'V',
          'ILE':'I','LEU':'L','MET':'M','PHE':'F','TYR':'Y',
          'TRP':'W'}

a2aa={'R':'ARG','H':'HIS','K':'LYS','D':'ASP','E':'GLU',
      'S':'SER','T':'THR','N':'ASN','Q':'GLN','C':'CYS',
      'U':'SEC','G':'GLY','P':'PRO','A':'ALA','V':'VAL',
      'I':'ILE','L':'LEU','M':'MET','F':'PHE','Y':'TYR',
      'W':'TRP'}

a2t={'D':'C','E':'C','H':'C','K':'C','R':'C',
         'P':'P','V':'H','M':'H','I':'H','L':'H',
         'F':'H','W':'H','G':'G','A':'H','C':'C',
         'T':'P','Q':'P','N':'P','Y':'P','S':'P'}

smaa1 =   ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X', '*']
smaa2 =   ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
aablast = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'J', 'Z', 'X', '*']
aaemboss = smaa1

b62edge = np.array([
        [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
        0.,  0.,  0.,  0.,  0.,  0.,  0., -2., -1.,  0., -4.],
        [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
        0.,  0.,  0.,  0.,  0.,  0.,  0., -1.,  0., -1., -4.],
        [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
        0.,  0.,  0.,  0.,  0.,  0.,  0.,  3.,  0., -1., -4.],
        [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
        0.,  0.,  0.,  0.,  0.,  0.,  0.,  4.,  1., -1., -4.],
        [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
        0.,  0.,  0.,  0.,  0.,  0.,  0., -3., -3., -2., -4.],
        [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
        0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  3., -1., -4.],
        [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
        0.,  0.,  0.,  0.,  0.,  0.,  0.,  1.,  4., -1., -4.],
        [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
        0.,  0.,  0.,  0.,  0.,  0.,  0., -1., -2., -1., -4.],
        [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
        0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0., -1., -4.],
        [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
        0.,  0.,  0.,  0.,  0.,  0.,  0., -3., -3., -1., -4.],
        [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
        0.,  0.,  0.,  0.,  0.,  0.,  0., -4., -3., -1., -4.],
        [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
        0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  1., -1., -4.],
        [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
        0.,  0.,  0.,  0.,  0.,  0.,  0., -3., -1., -1., -4.],
        [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
        0.,  0.,  0.,  0.,  0.,  0.,  0., -3., -3., -1., -4.],
        [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
        0.,  0.,  0.,  0.,  0.,  0.,  0., -2., -1., -2., -4.],
        [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
        0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0., -4.],
        [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
        0.,  0.,  0.,  0.,  0.,  0.,  0., -1., -1.,  0., -4.],
        [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
        0.,  0.,  0.,  0.,  0.,  0.,  0., -4., -3., -2., -4.],
        [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
        0.,  0.,  0.,  0.,  0.,  0.,  0., -3., -2., -1., -4.],
        [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
        0.,  0.,  0.,  0.,  0.,  0.,  0., -3., -2., -1., -4.],
        [-2., -1.,  3.,  4., -3.,  0.,  1., -1.,  0., -3., -4.,  0., -3.,
         -3., -2.,  0., -1., -4., -3., -3.,  4.,  1., -1., -4.],
        [-1.,  0.,  0.,  1., -3.,  3.,  4., -2.,  0., -3., -3.,  1., -1.,
         -3., -1.,  0., -1., -3., -2., -2.,  1.,  4., -1., -4.],
        [ 0., -1., -1., -1., -2., -1., -1., -1., -1., -1., -1., -1., -1.,
         -1., -2.,  0.,  0., -2., -1., -1., -1., -1., -1., -4.],
        [-4., -4., -4., -4., -4., -4., -4., -4., -4., -4., -4., -4., -4.,
         -4., -4., -4., -4., -4., -4., -4., -4., -4., -4.,  1.]])

b80blast = np.array([
        [ 5, -2, -2, -2, -1, -1, -1,  0, -2, -2, -2, -1, -1, -3, -1,  1,  0, -3, -2,  0, -2, -2, -1, -1, -6],
        [-2,  6, -1, -2, -4,  1, -1, -3,  0, -3, -3,  2, -2, -4, -2, -1, -1, -4, -3, -3, -1, -3,  0, -1, -6],
        [-2, -1,  6,  1, -3,  0, -1, -1,  0, -4, -4,  0, -3, -4, -3,  0,  0, -4, -3, -4,  5, -4,  0, -1, -6],
        [-2, -2,  1,  6, -4, -1,  1, -2, -2, -4, -5, -1, -4, -4, -2, -1, -1, -6, -4, -4,  5, -5,  1, -1, -6],
        [-1, -4, -3, -4,  9, -4, -5, -4, -4, -2, -2, -4, -2, -3, -4, -2, -1, -3, -3, -1, -4, -2, -4, -1, -6],
        [-1,  1,  0, -1, -4,  6,  2, -2,  1, -3, -3,  1,  0, -4, -2,  0, -1, -3, -2, -3,  0, -3,  4, -1, -6],
        [-1, -1, -1,  1, -5,  2,  6, -3,  0, -4, -4,  1, -2, -4, -2,  0, -1, -4, -3, -3,  1, -4,  5, -1, -6],
        [ 0, -3, -1, -2, -4, -2, -3,  6, -3, -5, -4, -2, -4, -4, -3, -1, -2, -4, -4, -4, -1, -5, -3, -1, -6],
        [-2,  0,  0, -2, -4,  1,  0, -3,  8, -4, -3, -1, -2, -2, -3, -1, -2, -3,  2, -4, -1, -4,  0, -1, -6],
        [-2, -3, -4, -4, -2, -3, -4, -5, -4,  5,  1, -3,  1, -1, -4, -3, -1, -3, -2,  3, -4,  3, -4, -1, -6],
        [-2, -3, -4, -5, -2, -3, -4, -4, -3,  1,  4, -3,  2,  0, -3, -3, -2, -2, -2,  1, -4,  3, -3, -1, -6],
        [-1,  2,  0, -1, -4,  1,  1, -2, -1, -3, -3,  5, -2, -4, -1, -1, -1, -4, -3, -3, -1, -3,  1, -1, -6],
        [-1, -2, -3, -4, -2,  0, -2, -4, -2,  1,  2, -2,  6,  0, -3, -2, -1, -2, -2,  1, -3,  2, -1, -1, -6],
        [-3, -4, -4, -4, -3, -4, -4, -4, -2, -1,  0, -4,  0,  6, -4, -3, -2,  0,  3, -1, -4,  0, -4, -1, -6],
        [-1, -2, -3, -2, -4, -2, -2, -3, -3, -4, -3, -1, -3, -4,  8, -1, -2, -5, -4, -3, -2, -4, -2, -1, -6],
        [ 1, -1,  0, -1, -2,  0,  0, -1, -1, -3, -3, -1, -2, -3, -1,  5,  1, -4, -2, -2,  0, -3,  0, -1, -6],
        [ 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -2, -1, -1, -2, -2,  1,  5, -4, -2,  0, -1, -1, -1, -1, -6],
        [-3, -4, -4, -6, -3, -3, -4, -4, -3, -3, -2, -4, -2,  0, -5, -4, -4, 11,  2, -3, -5, -3, -3, -1, -6],
        [-2, -3, -3, -4, -3, -2, -3, -4,  2, -2, -2, -3, -2,  3, -4, -2, -2,  2,  7, -2, -3, -2, -3, -1, -6],
        [ 0, -3, -4, -4, -1, -3, -3, -4, -4,  3,  1, -3,  1, -1, -3, -2,  0, -3, -2,  4, -4,  2, -3, -1, -6],
        [-2, -1,  5,  5, -4,  0,  1, -1, -1, -4, -4, -1, -3, -4, -2,  0, -1, -5, -3, -4,  5, -4,  0, -1, -6],
        [-2, -3, -4, -5, -2, -3, -4, -5, -4,  3,  3, -3,  2,  0, -4, -3, -1, -3, -2,  2, -4,  3, -3, -1, -6],
        [-1,  0,  0,  1, -4,  4,  5, -3,  0, -4, -3,  1, -1, -4, -2,  0, -1, -3, -3, -3,  0, -3,  5, -1, -6],
        [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -6],
        [-6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6, -6,  1]])

b62blast = np.array([
        [ 4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -2, -1, -1, -1, -4],
        [-1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -1, -2,  0, -1, -4],
        [-2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,  4, -3,  0, -1, -4],
        [-2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,  4, -3,  1, -1, -4],
        [ 0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -1, -3, -1, -4],
        [-1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,  0, -2,  4, -1, -4],
        [-1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,  1, -3,  4, -1, -4],
        [ 0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -1, -4, -2, -1, -4],
        [-2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,  0, -3,  0, -1, -4],
        [-1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3, -3,  3, -3, -1, -4],
        [-1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -4,  3, -3, -1, -4],
        [-1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2,  0, -3,  1, -1, -4],
        [-1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -3,  2, -1, -1, -4],
        [-2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -3,  0, -3, -1, -4],
        [-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -2, -3, -1, -1, -4],
        [ 1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,  0, -2,  0, -1, -4],
        [ 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -1, -1, -1, -1, -4],
        [-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -4, -2, -2, -1, -4],
        [-2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -3, -1, -2, -1, -4],
        [ 0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -3,  2, -2, -1, -4],
        [-2, -1,  4,  4, -3,  0,  1, -1,  0, -3, -4,  0, -3, -3, -2,  0, -1, -4, -3, -3,  4, -3,  0, -1, -4],
        [-1, -2, -3, -3, -1, -2, -3, -4, -3,  3,  3, -3,  2,  0, -3, -2, -1, -2, -1,  2, -3,  3, -3, -1, -4],
        [-1,  0,  0,  1, -3,  4,  4, -2,  0, -3, -3,  1, -1, -3, -1,  0, -1, -2, -2, -2,  0, -3,  4, -1, -4],
        [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -4],
        [-4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,  1]])

#  0        1                  2                     3           4                      5
# fullname, average resi mass, van der waals volume, polarities, number of heavy atoms  average accessible area
aadef={  #0        1       2  3  4           5
        'A':('Alanine',71.0788,67,9, 5,          31.5),
        'R':('Arginine',156.1876,148,15, 11,     93.8),
        'N':('Asparagine',114.1039,96,16, 8,     62.2),
        'D':('Aspartic.acid',115.0886,91,19, 8,  60.9),
        'C':('Cysteine',103.1448,86,7, 8,        13.9),
        'Q':('Glutamine',128.1308,114,17, 9,     74.0),
        'E':('Glutamic.acid',129.1155,109,18, 9, 72.3),
        'G':('Glycine',57.0520,48,11, 4,         25.2),
        'H':('Histidine',137.1412,118,10, 10,    46.7),
        'I':('Isoleucine',113.1595,124,1, 8,     23.0),
        'L':('Leucine',113.1595,124,3, 8,        29.0),
        'K':('Lysine',128.1742,135,20, 9,        110.3),
        'M':('Methionine',131.1986,124,5, 8,     30.5),
        'F':('Phenylalanine',147.1766,135,2, 11, 28.7),
        'P':('Proline',97.1167,90,13, 7,         53.7),
        'S':('Serine',87.0782,73,14, 6,          44.2),
        'T':('Threonine',101.1051,93,12, 7,      46.0),
        'W':('Tryptophan',186.2133,163,6, 14,    41.7),
        'Y':('Tyrosine',163.1760,141,8, 12,      59.1),
        'V':('Valine',99.1326,105,4, 7,          23.5),
        'U':('Cysteine',103.1448,86,7, 8,        13.9)
        }

'''
Column Headings:

1)  One-letter code (key)
2)  Three-letter code
3)  Name
4)  Chou-Fasman code for helix propensity
5)  Chou-Fasman code for sheet propensity
6)  Chou-Fasman helix propensity values
7)  Chou-Fasman sheet propensity values
8)  Amino acid molecular weight
9)  pKa value for free amino acid carboxylate
10) pKa value for free amino acid amine
11) pKa value for amino acid side chain
12) Number of carbon atoms in amino acid
13) Number of hydrogen atoms in amino acid zwitterion
14) Number of nitrogen atoms in amino acid
15) Number of oxygen atoms in amino acid
16) Number of sulfur atoms in amino acid
17) Area in the standard state (standard state accessibility is defined
        as the average surface area that residue has in a ensemble
        of Gly-X-Gly tripeptides)
18) Average accessible area in proteins
19) Average are buried upon transfer from the standard state to the folded protein.
20) Mean fractional area loss, equal to the average area buried normalized
        by the standard state area
21) Residue mass
22) Monoisotopic Mass
23) Total number of heavy atom
'''
aaprop = {
        'A':('ALA' ,  'Alanine'       ,  'H',  'I'   ,1.45  ,0.97   ,89.09  ,2.3 ,9.9   ,0    ,3 ,7 ,1,2,0,118.1,31.5 ,86.6 ,.74,71.08 ,71.03711 , 6),
        'C':('CYS' ,  'Cysteine'      ,  'i',  'h'   ,0.77  ,1.30   ,121.16 ,1.8 ,10.8  ,8.65 ,3 ,7 ,1,2,1,146.1,13.9 ,132.3,.91,103.14,103.00919, 7),
        'D':('ASP' ,  'Aspartic Acid' ,  'i',  'i'   ,0.98  ,0.80   ,133.10 ,2.0 ,10.0  ,4.04 ,4 ,7 ,1,4,0,158.7,60.9 ,97.8 ,.62,115.09,115.02694, 9),
        'E':('GLU' ,  'Glutamic Acid' ,  'H',  'B'   ,1.53  ,0.26   ,147.13 ,2.2 ,9.7   ,4.39 ,5 ,9 ,1,4,0,186.2,72.3 ,113.9,.62,129.12,129.04259,10),
        'F':('PHE' ,  'Phenylalanine' ,  'h',  'h'   ,1.12  ,1.28   ,165.19 ,1.8 ,9.1   ,0    ,9 ,11,1,2,0,222.8,28.7 ,194.1,.88,147.18,147.06841,12),
        'G':('GLY' ,  'Glycine'       ,  'B',  'i'   ,0.53  ,0.81   ,75.07  ,2.4 ,9.8   ,0    ,2 ,5 ,1,2,0,88.1 ,25.2 ,62.9 ,.72,57.05 ,57.02146 , 5),
        'H':('HIS' ,  'Histidine'     ,  'h',  'b'   ,1.24  ,0.71   ,155.16 ,1.8 ,9.2   ,6.75 ,6 ,9 ,3,2,0,202.5,46.7 ,155.8,.78,137.14,137.05891,11),
        'I':('ILE' ,  'Isoleucine'    ,  'I',  'H'   ,1.00  ,1.60   ,131.17 ,2.4 ,9.7   ,0    ,6 ,13,1,2,0,181  ,23   ,158  ,.88,113.16,113.08406, 9),
        'K':('LYS' ,  'Lysine'        ,  'I',  'b'   ,1.07  ,0.74   ,146.19 ,2.2 ,9.2   ,11.0 ,6 ,14,2,2,0,225.8,110.3,115.5,.52,128.17,128.09496,10),
        'L':('LEU' ,  'Leucine'       ,  'H',  'h'   ,1.34  ,1.22   ,131.17 ,2.4 ,9.60  ,0    ,6 ,13,1,2,0,193.1,29   ,164.1,.85,113.16,113.08406, 9),
        'M':('MET' ,  'Methionine'    ,  'h',  'H'   ,1.20  ,1.67   ,149.21 ,2.3 ,9.2   ,0    ,5 ,11,1,2,1,203.4,30.5 ,172.9,.85,131.19,131.04049, 9),
        'N':('ASN' ,  'Asparagine'    ,  'b',  'b'   ,0.73  ,0.65   ,132.12 ,2.0 ,8.8   ,0    ,4 ,8 ,2,3,0,165.5,62.2 ,103.3,.63,114.10,114.04293, 9),
        'P':('PRO' ,  'Proline'       ,  'B',  'b'   ,0.59  ,0.62   ,115.13 ,2.0 ,10.6  ,0    ,5 ,9 ,1,2,0,146.8,53.7 ,92.9 ,.64,97.12 ,97.05276 , 8),
        'Q':('GLN' ,  'Glutamine'     ,  'h',  'h'   ,1.17  ,1.23   ,146.15 ,2.2 ,9.1   ,0    ,5 ,10,2,3,0,193.2,74   ,119.2,.62,128.13,128.05858,10),
        'R':('ARG' ,  'Arginine'      ,  'i',  'i'   ,0.79  ,0.90   ,174.20 ,1.8 ,9.0   ,12.5 ,6 ,14,4,2,0,256  ,93.8 ,162.2,.64,156.19,156.10111,12),
        'S':('SER' ,  'Serine'        ,  'i',  'b'   ,0.79  ,0.72   ,105.09 ,2.1 ,9.2   ,0    ,3 ,7 ,1,3,0,129.8,44.2 ,85.6 ,.66,87.08 ,87.03203 , 7),
        'T':('THR' ,  'Threonine'     ,  'i',  'h'   ,0.82  ,1.20   ,119.12 ,2.6 ,10.4  ,0    ,4 ,9 ,1,3,0,152.5,46   ,106.5,.70,101.11,101.04768, 8),
        'V':('VAL' ,  'Valine'        ,  'h',  'H'   ,1.14  ,1.65   ,117.15 ,2.3 ,9.6   ,0    ,5 ,11,1,2,0,164.5,23.5 ,141  ,.86,99.13 ,99.06841 , 8),
        'W':('TRP' ,  'Tryptophan'    ,  'h',  'h'   ,1.14  ,1.19   ,204.22 ,2.4 ,9.4   ,0    ,11,12,2,2,0,266.3,41.7 ,224.6,.85,186.21,186.07931,15),
        'Y':('TYR' ,  'Tyrosine'      ,  'b',  'h'   ,0.61  ,1.29   ,181.19 ,2.20,9.1   ,9.75 ,9 ,11,1,3,0,236.8,59.1 ,177.7,.76,163.18,163.06333,13),
        '.':('gap' ,  'gap'           ,  '.' , '.'   ,0     ,0      ,0      ,0   ,0     ,0    ,0 ,0 ,0,0,0,0    ,0    ,0    ,0  ,0     ,0        ,0 ),
        '-':('gap' ,  'gap'           ,  '.' , '.'   ,0     ,0      ,0      ,0   ,0     ,0    ,0 ,0 ,0,0,0,0    ,0    ,0    ,0  ,0     ,0        ,0 ),
        'O':('gap' ,  'gap'           ,  '.' , '.'   ,0     ,0      ,0      ,0   ,0     ,0    ,0 ,0 ,0,0,0,0    ,0    ,0    ,0  ,0     ,0        ,0 ),
        'Z':('gap' ,  'gap'           ,  '.' , '.'   ,0     ,0      ,0      ,0   ,0     ,0    ,0 ,0 ,0,0,0,0    ,0    ,0    ,0  ,0     ,0        ,0 ),
        'X':('gap' ,  'gap'           ,  '.' , '.'   ,0     ,0      ,0      ,0   ,0     ,0    ,0 ,0 ,0,0,0,0    ,0    ,0    ,0  ,0     ,0        ,0 ),
        'U':('gap' ,  'gap'           ,  '.' , '.'   ,0     ,0      ,0      ,0   ,0     ,0    ,0 ,0 ,0,0,0,0    ,0    ,0    ,0  ,0     ,0        ,0 ),
        'J':('gap' ,  'gap'           ,  '.' , '.'   ,0     ,0      ,0      ,0   ,0     ,0    ,0 ,0 ,0,0,0,0    ,0    ,0    ,0  ,0     ,0        ,0 ),
        'B':('gap' ,  'gap'           ,  '.' , '.'   ,0     ,0      ,0      ,0   ,0     ,0    ,0 ,0 ,0,0,0,0    ,0    ,0    ,0  ,0     ,0        ,0 )
}

# MSA score
aapscore = {
        'ssp': {'.':0, 'H':1, 'I':2, 'B':3, 'h':4, 'i':5, 'b':6},
        'hat': {0:0, 5:1, 6:2, 7:3, 8:4, 9:5, 10:6, 11:7, 12:8, 13:9, 15:10}
}

# kidera factors
# mean: 0, std: 1
# gap: 3
kidera = {
        'A':(-1.56, -1.67, -0.97, -0.27, -0.93, -0.78,  -0.2, -0.08,  0.21, -0.48),
        'R':( 0.22,  1.27,  1.37,  1.87,  -1.7,  0.46,  0.92, -0.39,  0.23,  0.93),
        'N':( 1.14, -0.07, -0.12,  0.81,  0.18,  0.37, -0.09,  1.23,   1.1, -1.73),
        'D':( 0.58, -0.22, -1.58,  0.81, -0.92,  0.15, -1.52,  0.47,  0.76,   0.7),
        'C':( 0.12, -0.89,  0.45, -1.05, -0.71,  2.41,  1.52, -0.69,  1.13,   1.1),
        'Q':(-0.47,  0.24,  0.07,   1.1,   1.1,  0.59,  0.84, -0.71, -0.03, -2.33),
        'E':(-1.45,  0.19, -1.61,  1.17, -1.31,   0.4,  0.04,  0.38, -0.35, -0.12),
        'G':( 1.46, -1.96, -0.23, -0.16,   0.1, -0.11,  1.32,  2.36, -1.66,  0.46),
        'H':(-0.41,  0.52, -0.28,  0.28,  1.61,  1.01, -1.85,  0.47,  1.13,  1.63),
        'I':(-0.73, -0.16,  1.79, -0.77, -0.54,  0.03, -0.83,  0.51,  0.66, -1.78),
        'L':(-1.04,    0,  -0.24,  -1.1, -0.55, -2.05,  0.96, -0.76,  0.45,  0.93),
        'K':(-0.34,  0.82, -0.23,   1.7,  1.54, -1.62,  1.15, -0.08, -0.48,   0.6),
        'M':(-1.4,   0.18, -0.42, -0.73,     2,  1.52,  0.26,  0.11, -1.27,  0.27),
        'F':(-0.21,  0.98, -0.36, -1.43,  0.22, -0.81,  0.67,   1.1,  1.71, -0.44),
        'P':( 2.06, -0.33, -1.15, -0.75,  0.88, -0.45,   0.3,  -2.3,  0.74, -0.28),
        'S':( 0.81, -1.08,  0.16,  0.42, -0.21, -0.43, -1.89, -1.15, -0.97, -0.23),
        'T':( 0.26,  -0.7,  1.21,  0.63,  -0.1,  0.21,  0.24, -1.15, -0.56,  0.19),
        'W':(  0.3,   2.1, -0.72, -1.57, -1.16,  0.57, -0.48,  -0.4,  -2.3,  -0.6),
        'Y':( 1.38,  1.48,   0.8, -0.56,     0, -0.68, -0.31,  1.03, -0.05,  0.53),
        'V':(-0.74, -0.71,  2.04,  -0.4,   0.5, -0.81, -1.07,  0.06, -0.46,  0.65),
        '.':( 3.00,  3.00,  3.00,  3.00,  3.00,  3.00,  3.00,  3.00,  3.00,  3.00),
        '-':( 3.00,  3.00,  3.00,  3.00,  3.00,  3.00,  3.00,  3.00,  3.00,  3.00),
        'O':( 3.00,  3.00,  3.00,  3.00,  3.00,  3.00,  3.00,  3.00,  3.00,  3.00),
        'Z':( 3.00,  3.00,  3.00,  3.00,  3.00,  3.00,  3.00,  3.00,  3.00,  3.00),
        'X':( 3.00,  3.00,  3.00,  3.00,  3.00,  3.00,  3.00,  3.00,  3.00,  3.00),
        'U':( 3.00,  3.00,  3.00,  3.00,  3.00,  3.00,  3.00,  3.00,  3.00,  3.00),
        'J':( 3.00,  3.00,  3.00,  3.00,  3.00,  3.00,  3.00,  3.00,  3.00,  3.00),
        'B':( 3.00,  3.00,  3.00,  3.00,  3.00,  3.00,  3.00,  3.00,  3.00,  3.00)
}


# print repr(dict((k,cp.aapscore['ssp'][cp.aaprop[k][2]]) for k in cp.aaprop))
aascore = {
        'aa' :  {
                        'J':0, 'O':0, 'Z':0, 'U':0,'X':0,'-': 0,'.': 0,'A': 1,'C': 2,'D': 3,'E': 4,'F': 5,'G': 6,'H': 7,'I': 8,'K': 9,
                        'L': 10,'M': 11,'N': 12,'P': 13,'Q': 14,'R': 15,'S': 16,'T': 17,'V': 18,'W': 19,'Y': 20, 'B': 0
                        },

        'ssp' : {
                        '-': 0, '.': 0, 'A': 1, 'C': 5, 'E': 1, 'D': 5, 'G': 3, 'F': 4, 'I': 2, 'H': 4, 'K': 2, 'J': 0, 'M': 4, 'B': 0,
                        'L': 1, 'O': 0, 'N': 6, 'Q': 4, 'P': 3, 'S': 5, 'R': 5, 'U': 0, 'T': 5, 'W': 4, 'V': 4, 'Y': 6, 'X': 0, 'Z': 0
                        },

        'hat' : {
                        '-': 0, '.': 0, 'A': 2, 'C': 3, 'B': 0, 'E': 6, 'D': 5, 'G': 1, 'F': 8, 'I': 5, 'H': 7, 'K': 6, 'J': 0, 'M': 5,
                        'L': 5, 'O': 0, 'N': 5, 'Q': 6, 'P': 4, 'S': 3, 'R': 8, 'U': 0, 'T': 4, 'W': 10, 'V': 4, 'Y': 9, 'X': 0, 'Z': 0
                        },

        'gthat':{ # for fdr
                        'B':0, 'J':0, 'O':0, 'Z':0, 'U':0, 'X':0, '-':0, '.':0,
                        'A':3, 'G':3,
                        'C':2, 'S':2, 'T':2, 'N':2, 'Q':2,
                        'D':1, 'R':1, 'E':1, 'K':1, 'H':1,
                        'F':5, 'Y':5, 'W':5,
                        'P':4, 'V':4, 'I':4, 'L':4, 'M':4
                        },
        'bin':  {
                        'J':0, 'O':0, 'Z':0, 'U':0,'X':0,'-': 0,'.': 0,'A': 1,'C': 1,'D': 1,'E': 1,'F': 1,'G': 1,'H': 1,'I': 1,'K': 1,
                        'L': 1,'M': 1,'N': 1,'P': 1,'Q': 1,'R': 1,'S': 1,'T': 1,'V': 1,'W': 1,'Y': 1, 'B': 0
                        },
        'triple':       {
                        'J':0, 'O':0, 'Z':0, 'U':0,'X':0,'-': 0,'.': 0.5,'A': 1,'C': 1,'D': 1,'E': 1,'F': 1,'G': 1,'H': 1,'I': 1,'K': 1,
                        'L': 1,'M': 1,'N': 1,'P': 1,'Q': 1,'R': 1,'S': 1,'T': 1,'V': 1,'W': 1,'Y': 1, 'B': 0
                        },
        # from "Simplified amino acid alphabets for protein fold recognition and implications for folding" 2000
        'r2000_2' :     {
                        'J':0, 'O':0, 'Z':0, 'U':0,'X':0,'-': 0,'.': 0,'A': 1,'C': 1,'D': 2,'E': 2,'F': 1,'G': 1,'H': 2,'I': 1,'K': 2,
                        'L': 1,'M': 1,'N': 2,'P': 1,'Q': 2,'R': 2,'S': 1,'T': 1,'V': 1,'W': 1,'Y': 1, 'B': 0
                        },
        'r2000_4' :     {
                        'J':0, 'O':0, 'Z':0, 'U':0,'X':0,'-': 0,'.': 0,'A': 2,'C': 1,'D': 4,'E': 4,'F': 3,'G': 2,'H': 4,'I': 1,'K': 4,
                        'L': 1,'M': 1,'N': 4,'P': 2,'Q': 4,'R': 4,'S': 2,'T': 2,'V': 1,'W': 3,'Y': 3, 'B': 0
                        },
        'r2000_8' :     {
                        'J':0, 'O':0, 'Z':0, 'U':0,'X':0,'-': 0,'.': 0,'A': 2,'C': 1,'D': 6,'E': 6,'F': 5,'G': 2,'H': 8,'I': 1,'K': 7,
                        'L': 1,'M': 1,'N': 6,'P': 4,'Q': 6,'R': 7,'S': 3,'T': 3,'V': 1,'W': 5,'Y': 5, 'B': 0
                        },
        'r2000_10' :    {
                        'J':0, 'O':0, 'Z':0, 'U':0,'X':0,'-': 0,'.': 0,'A': 3,'C': 2,'D': 8,'E': 8,'F': 7,'G': 4,'H': 10,'I': 1,'K': 9,
                        'L': 1,'M': 1,'N': 8,'P': 6,'Q': 8,'R': 9,'S': 5,'T': 5,'V': 1,'W': 7,'Y': 7, 'B': 0
                        },
        'r2000_15' :    {
                        'J':0, 'O':0, 'Z':0, 'U':0,'X':0,'-': 0,'.': 0,'A': 3,'C': 2,'D': 11,'E': 10,'F': 8,'G': 4,'H': 15,'I': 1,'K': 14,
                        'L': 1,'M': 1,'N': 12,'P': 7,'Q': 13,'R': 14,'S': 5,'T': 6,'V': 1,'W': 9,'Y': 8, 'B': 0
                        },
        'rna' : {
                '-':0, '.':0, 'm':0, 'M':0, 'n':0, 'N':0, 'r':0, 'R':0, 'w':0, 'W':0, 'y':0, 'Y':0, 'k':0, 'K':0,
                's':0, 'S':0, 'v':0, 'V':0, 'h':0, 'H':0, 'd':0, 'D':0, 'b':0, 'B':0,
                'a':1, 'A':1, 'c':2, 'C':2, 'g':3, 'G':3, 'u':4, 'U':4, 't':5, 'T':5,
                },                        
        'null' : {}
}

# print repr(dict((cp.aascore['aa'][k], k) for k in cp.aascore['aa'] if cp.aascore['aa'][k]!=0))
# reverse dictionary of aascore
scoreaa = {
        'aa' : {
                        0: '.', 1: 'A', 2: 'C', 3: 'D', 4: 'E', 5: 'F', 6: 'G', 7: 'H', 8: 'I', 9: 'K', 10: 'L', 11: 'M', 12: 'N',
                        13: 'P', 14: 'Q', 15: 'R', 16: 'S', 17: 'T', 18: 'V', 19: 'W', 20: 'Y'
                        }
}

scorerdict = {
                                0:'.',1:'A',2:'C',3:'D',4:'E',5:'F',6:'G',7:'H',8:'I',9:'K',
                                10:'L',11:'M',12:'N',13:'P',14:'Q',15:'R',16:'S',17:'T',18:'V',19:'W',20:'Y'
                        }

msaaa = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M',
                'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X', '.', '-']

chargedaa = ['R', 'H', 'K', 'D', 'E']
aacharge = {'R':'P', 'H':'P', 'K':'P', 'D':'N', 'E':'N'}

# EBLOSUM matrices gaps preset
#gapb80 = [(9,2), (8,2), (7,2), (6,2), (11,1), (10,1), (9,1)]
gapb62=[(9,2), (8,2), (7,2), (6,2), (11,1), (10,1), (9,1), (11, 2), (10, 2), (13, 1), (12, 1)]
gapb45=[(13, 3),(12, 3),(11, 3),(10, 3),(16, 2),(15, 2),(14, 2),(13, 2),(12, 2),(19, 1),(18, 1),(17, 1),(16, 1)]
gapb50=[(13, 3),(12, 3),(11, 3),(10, 3),(9, 3),(16, 2),(15, 2),(14, 2),(13, 2),(12, 2),(19, 1),(18, 1),(17, 1),(16, 1),(15, 1)]
gapb80=[(25, 2),(13, 2),(9, 2),(8, 2),(7, 2),(6, 2),(11, 1),(10, 1),(9, 1)]
gapb90=[(9, 2),(8, 2),(7, 2),(6, 2),(11, 1),(10, 1),(9, 1)]
gappam30=[(7, 2),(6, 2),(5, 2),(10, 1),(9, 1),(8, 1),(15, 3),(14, 2),(14, 1),(13, 3)]
gappam70=[(8, 2),(7, 2),(6, 2),(11, 1),(10, 1),(9, 1),(11, 2),(12, 3)]
gappam250=[(15, 3),(14, 3),(13, 3),(12, 3),(11, 3),(17, 2),(16, 2),(15, 2),(14, 2),(13, 2),(21, 1),(20, 1),(19, 1),(18, 1),(17, 1)]

gapdict = {
        'BLOSUM80':gapb80,
        'BLOSUM62':gapb62,
        'BLOSUM50':gapb50,
        'BLOSUM45':gapb45,
        'PAM250':gappam250,
        'BLOSUM90':gapb90,
        'PAM30':gappam30,
        'PAM70':gappam70
}

colormap8 = ['[0.7773,0.1328,0.1484]',
                        '[0.3125,0.6055,0.7188]',
                        '[0.9727,0.8086,0.7461]',
                        '[0.2016,0.2055,0.2211]',
                        '[0.9805,0.5469,0.1953]',
                        '[0.3750,0.8984,0.7539]',
                        '[0.9242,0.7015,0.3305]',
                        '[0.0391,0.5508,0.0430]']

time0 = time.time()

def _fatal():
    exit(1)

def _warning():
    return False
# exit causes problem in mprun
def _err(msg, errcallback=_fatal):
    curframe = inspect.currentframe()
    calframe = inspect.getouterframes(curframe, 1)
    print('[err:%s:%s()] %s' % (' '.join(sys.argv),calframe[1][3], msg))
    errcallback()

def _info(msg, flag='INFO'):
    global time0
    time1 = time.time()
    #tick = (time1-time0)/3600
    tick = int((time1-time0))
    curframe = inspect.currentframe()
    calframe = inspect.getouterframes(curframe, 1)
    # time pid tick flag/level msg
    info = '%s|%d|%d|%s|%s' % (time.strftime('%Y-%m-%d %H:%M:%S'), os.getpid(), tick, flag, msg)
    #info = '%d:%s:%s()' % (os.getpid(),':'.join(sys.argv[1:]), calframe[1][3])
    time0 = time1
    print(info)





# return ith column of matrix (list(list))
#       [...]
# [ [...] ]
#       [...]
def column(mat, i):
    return [row[i] for row in mat]

# calculate the Euclidean distance between two vectors
def dist(v1, v2):
    return np.linalg.norm(np.array(v1)-np.array(v2))


# mp constant definition
mp_info = 0
mp_log = 1
mp_checkin = 2

def dispatch(module):
    return getattr(sys.modules[module], sys.argv[1])(sys.argv[2:]) if (len(sys.argv) >= 2 and sys.argv[1] in dir(sys.modules[module])) else _err('cmd not found')


# used in utils_mprun.py
def dcall(callstr):
    strarr = callstr.split()
    modu = strarr[1][:-3] # to ignore '.py'
    func = strarr[2]
    param = strarr[3:]

    ins_func = getattr(__import__(modu), func)
    return ins_func(param)

# used in utils_mprun.py
def drun(runstr):
    return subprocess.Popen(runstr, stdout=subprocess.PIPE, shell=True).communicate()[0].strip()

# calculate adjusted mean and std
# introduced in sdii shadow paper
def adjmean(nlist, c):
	nplist = np.array(nlist)
	if c==0: 
		return nplist.mean()
	outlier = nplist.mean()+c*nplist.std()
	nplist_no_outlier = np.array([v for v in nlist if v < outlier])
	adjm = nplist_no_outlier.mean()
	return adjm    

# transfer a skewed distribution into normal
# input: vec, a list of float values
# input: method: log, rsquare, csquare, default = log
def disttrans(vec, method='log'):
    if method == 'csquare':
        op = lambda x: x ** ( 1. / 3 )
    elif method == 'rsquare':
        op = lambda x: x ** ( 1. / 2 )
    else:
        print('use default log transformation')
        op = lambda x: math.log(x,2)
    return [op(v) for v in vec]

# set /0 = 0
def div0(a, b):
    from numpy import errstate, true_divide, isfinite, isscalar
    with errstate(divide='ignore', invalid='ignore'):
        c = true_divide(a, b)
        if isscalar(c):
            if not isfinite(c):
                c = 0
        else:
            c[~isfinite(c)] = 0.  # -inf inf NaN
    return c

# (joint) entropy calculation
# input: a list of np.array()
def entropy(X):
    return np.sum(-p * np.log2(p) if p > 0 else 0 for p in (np.mean(reduce(np.logical_and, (predictions == c for predictions, c in zip(X, classes)))) for classes in itertools.product(*[set(x) for x in X])))


# generate fasta entries from a fasta file
def fasta_iter(fastafile):
    fh = open(fastafile)
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        header = next(header)[1:].strip()
        seq = "".join(s.strip() for s in next(faiter))
        yield header, seq


# return the freuency table for a given string
def freq(word):
    # 2.7+
    # return letters = collections.Counter('google')
    letters = collections.defaultdict(int)
    for letter in word:
        letters[letter] += 1
    return letters


# X: data[:,colset].T
#   0   1   2   3   4
#[[ 1.  1.  3.  2.  4.]
# [ 1.  1.  3.  1.  1.]]
# output:
# (1.0, 1.0), [0, 1]
# (2.0, 1.0), [3]
# (3.0, 3.0), [2]
# (4.0, 1.0), [4]
def freqlookup(X):
    lookup = []
    # generate cartesian product from a list (* is for a list of lists)
    for classes in itertools.product(*[set(x) for x in X]):
        #print 'start iteration:'
        #print ('c in [set([1.0, 2.0, 3.0, 4.0]), set([1.0, 3.0])] : ', classes)
        #print 'zip x and c: ', zip(X, classes)
        #print repr([predictions == c for predictions, c in zip(X, classes)])
        f = reduce(np.logical_and, (predictions == c for predictions, c in zip(X, classes)))
        #print 'normal freq for %s : %.02f' % (repr(classes), sum(f))
        #print 'weight freq for %s : %.02f' % (repr(classes), sum(w[[i for i in xrange(len(f)) if f[i]!=False]]))
        idx = [i for i in range(len(f)) if f[i]!=False]
        if len(idx)!=0:
            lookup.append((classes, idx))
    return lookup



def getPDBUniprotMap(mapfile):
    # duplicate from utils_msa.py
    # called in utils_msa.py ncg2sdiicol()
    # load map between pdb residue ID and MSA uniprot position ID
    # dictionary element: ('A9', (14, 'V')) : (chain+resi, (MSA position index, resn))
    # mapfile:
    # AT284 1218 T  : chain A residue T resn 284 => position 1218 (start from 0) resn T
    # AE285 1220 e  : lowercase exists!
    # AR286 -1 -    : residue number that cannot map to MSA position (does not exist)
    posmap = {}
    with open(mapfile) as fp:
        for line in fp:
            line = line.strip()
            if len(line) < 1:
                print('getPDBUniprotMap: error loading map: %s' % mapfile)
            strArray = line.split(' ')
            key = strArray[0][0] + strArray[0][2:]
            value = (int(strArray[1]), strArray[2].upper())
            posmap[key] = value
    print('getPDBUniprotMap: %s %d maps loaded' % (mapfile, len(posmap)))
    return posmap

# return the percentage of how many identical characters
def hamming_similarity(s1,s2):
    assert len(s1) == len(s2), 'length not match: \n%s\n%s' % (s1,s2)
    return 100*sum(ch1 == ch2 for ch1, ch2 in zip(s1, s2))/float(len(s1))

def hamming_diff(s1,s2):
    assert len(s1) == len(s2), 'length not match: \n%s\n%s' % (s1,s2)
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

# x: str list in np.array format
#$ python proc_hammingweight.py t_hamming_weight.score 0.9
#array([1, 2, 2, 3, 4, 1, 1], dtype=int32)
#defaultdict(<type 'int'>, {1: 3, 2: 2, 3: 1, 4: 1})
# 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 : 0.333
# 1, 1, 1, 1, 1, 1, 1, 1, 3, 3 : 0.500
# 1, 1, 1, 1, 1, 1, 1, 1, 3, 3 : 0.500
# 1, 1, 1, 1, 1, 3, 3, 3, 4, 4 : 1.000
# 1, 1, 1, 1, 1, 3, 3, 3,15,15 : 1.000
# 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 : 0.333
# 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 : 0.333
# max_d: cutoff distance, upper limit of the maximum ditance of two elemetns in the same cluster
#                                                  idx1       idx2        dist       # of elements
# linkage matrix Z elements: [ 14.     ,  79.     ,   0.05914,   2.     ]
# All indices idx >= len(X) actually refer to the cluster formed in Z[idx - len(X)]
'''
def hamming_weight(x, max_d):
    linkage_matrix = linkage(x, "single", metric='hamming')
    clusters = fcluster(linkage_matrix, max_d, criterion='distance')
    normdict = freq(clusters)
    return [(1.0/normdict[k]) for k in clusters]

# return a list of len(rows of x), clusters[i] = 'cluster ID which i belongs'
def hamming_cluster(x, max_d):
    linkage_matrix = linkage(x, "single", metric='hamming')
    clusters = fcluster(linkage_matrix, max_d, criterion='distance')
    return clusters
'''
# convert square formed matrix to flat list form
# names is a list of column names ordered according to matrix column indices
# output: list of vec3, c1 c2 value
def mat2flat(mat, names, sym=True):
    assert mat.shape[0] == mat.shape[1], 'cp.mat2flat::error:non-square matrix %s' % repr(mat.shape)
    return [(names[i],names[j],mat[i,j]) for i in range(mat.shape[0]) for j in range(i,mat.shape[0])] if sym==True else [(names[i],names[j],mat[i,j]) for i in range(mat.shape[0]) for j in range(mat.shape[0])]


# jaccard distance for two sets
def jaccard(a, b):
    c = a.intersection(b)
    print(repr(a), repr(b))
    print(repr(c))
    return 1 - (float(len(c)) / (len(a) + len(b) - len(c)))


# calculate value of n choose r
def ncr(n, r):
    return 1 if min(r, n-r) <= 0 else reduce(op.mul, range(n, n-r, -1))//reduce(op.mul, range(1, r+1))


# return a index set of n choose m
def ncrset(varnum, order):
    return [s for s in set(itertools.combinations(list(range(varnum)), order))]


def ncrvar(varset, order):
    return [s for s in set(itertools.combinations(varset, order))]

def normminmax(nplist):
    return 1.0*(nplist - np.min(nplist))/(np.max(nplist)-np.min(nplist))

def normmax(nplist):
    return 1.0*(nplist)/np.max(nplist)

# ordinary least square
# beta = (X^T*X)^(-1)*X^T*y
# input x: np.array (matrix)
# input y: np.array response
# return (m, b) # slope and intercept
def ols(X, y):
	return np.linalg.lstsq(np.vstack([X, np.ones(len(X))]).T, y, rcond=None)    

# given two strings
# normal sequence & aligned sequence
# return map 1. key=1  pos[s1] = s2; 2. key=2 pos[s2] = s1
# s1: aligned string index, s2: pdb sequence index
def posmap(s1, s2, key=1):
    gap = ['.', '-', '_']
    ps1 = s1.translate({ord(i):None for i in gap})
    ps2 = s2.translate({ord(i):None for i in gap})
    #print 'ps1: %s\nps2: %s' % (ps1, ps2)

    retmap={}
    if ps1!=ps2:
        print('error: not homo-str')
        print('ps1: %s\nps2: %s' % (ps1, ps2))
        return retmap

    i=0
    j=0
    while(i<len(s1) and j<len(s2)):
        if s1[i] in gap:
            i+=1
            continue
        if s2[j] in gap:
            j+=1
            continue
        if s1[i]==s2[j]:
            if key == 1:
                retmap[i] = j
            else:
                retmap[j] = i
            i+=1
            j+=1

    if len(retmap)!=len(ps1):
        print('error: incomplete map: len:%d, ps len: %d' % (len(retmap), len(ps1)))
        return False

    return retmap


# given two strings
# normal sequence & aligned sequence
# return map 1. key=1  pos[s1] = s2; 2. key=2 pos[s2] = s1
# s1: aligned string index, s2: pdb sequence index
def posmap1(s1, s2, key=1):
    gap = ['.', '-', '_']
    ps1 = s1.translate({ord(i):None for i in gap})
    ps2 = s2.translate({ord(i):None for i in gap})

    retmap={}
    if ps1!=ps2:
        print('error: not homo-str')
        print('ps1: %s\nps2: %s' % (ps1, ps2))
        return retmap
    i=0
    j=0
    while(i<len(s1) and j<len(s2)):
        if s1[i] in gap:
            i+=1
            continue
        if s2[j] in gap:
            j+=1
            continue
        if s1[i]==s2[j]:
            if key == 1:
                retmap[i] = j
            else:
                retmap[j] = i
            i+=1
            j+=1
    if len(retmap)!=len(ps1):
        print('error: incomplete map: len:%d, ps len: %d' % (len(retmap), len(ps1)))
        return False

    return retmap


# map between two gap extened sequences with the same original sequence
# normal sequence or aligned sequence
# return map pos[s1_index] = s2_index;
# index starts from 0
def posmap_homoseq(s1, s2):
    ps1 = s1.translate({ord(i):None for i in gap})
    ps2 = s2.translate({ord(i):None for i in gap})

    if ps1!=ps2:
        _err(_fatal, 'unmatched raw sequence\ns1: %s\ns2: %s\n' % (ps1, ps2))
    i=0
    j=0
    retmap=[]
    while(i<len(s1) and j<len(s2)):
        if s1[i] in gaps:
            i+=1
            continue
        if s2[j] in gaps:
            j+=1
            continue
        if s1[i]==s2[j]:
            retmap.append((i, j))
            i+=1
            j+=1

    if len(retmap)!=len(ps1):
        print('error: incomplete map: len:%d, ps len: %d' % (len(retmap), len(ps1)))
        return False
    return retmap


# extend version of posmap_homoseq()
# map between two gap extened sequences with one substring to the another without gaps
#       after pfamscan, a MSA sequence will be trimmed and insert gaps in different positions
# normal sequence or aligned sequence
# return map pos[s1_index] = s2_index;
# index starts from 0
def posmap_subseq(s1, s2):
    s1=s1.upper()
    s2=s2.upper()
    ps1 = s1.translate({ord(i):None for i in gaps})
    ps2 = s2.translate({ord(i):None for i in gaps})
    #print 'ps1:\n%s\n' % ps1
    #print 'ps2:\n%s\n' % ps2

    reverse, idx_set = (False, subseq_align(s1, s2, ps1.find(ps2))) if len(ps1) >= len(ps2) else (True, subseq_align(s2, s1, ps2.find(ps1)))

    if reverse == True:
        retmaplist = [(v, k) for (k, v) in idx_set]
    else:
        retmaplist = idx_set
    #retmap = dict((k,v) for (k,v) in retmaplist)
    return retmaplist

# update from posmap_subseq
# used in proc_psicov2msai.py
# return a dict instead of a list of tuples
def posmap_subseq_d(s1, s2):
    s1=s1.upper()
    s2=s2.upper()
    ps1 = s1.translate({ord(i):None for i in gaps})
    ps2 = s2.translate({ord(i):None for i in gaps})
    #print 'ps1:\n%s\n' % ps1
    #print 'ps2:\n%s\n' % ps2

    reverse, idx_set = (False, subseq_align(s1, s2, ps1.find(ps2))) if len(ps1) >= len(ps2) else (True, subseq_align(s2, s1, ps2.find(ps1)))

    if reverse == True:
        retmaplist = [(v, k) for (k, v) in idx_set]
    else:
        retmaplist = idx_set
    retmap = dict((k,v) for (k,v) in retmaplist)
    return retmap


# pairwise substitution
# return unified key
#def quad_permu(pair1, pair2):
# quad: list : ['A', 'C', 'D', 'G']
def quad_permu(quad):
    '''
    rank = 0  rank = 1  rank = 2  rank = 3
    A 0  C 1  C 0  A 1  D 0  G 1  G 0  D 1
    D 2  G 3  G 2  D 3  A 2  C 3  C 2  A 3
    '''
    rank = quad.index(min(quad))
    if rank == 1:
        quad[0],quad[1]=quad[1],quad[0]
        quad[2],quad[3]=quad[3],quad[2]
    elif rank == 2:
        quad[0],quad[2] = quad[2],quad[0]
        quad[1],quad[3] = quad[3],quad[1]
    elif rank == 3:
        quad[0],quad[3] = quad[3],quad[0]
        quad[2],quad[1] = quad[1],quad[2]
    return ''.join(quad)


# return which type of pair substitution the quadstr is
# quadstr = 'ACDG'
def quadtype(quadstr):
    count = 0
    pA = quadstr[0:2] # 'AC'
    pB = quadstr[2:]  # 'DG'
    #print (pA,pB)
    for i in [0,1]:
        if pA[i]!=pB[i]:
            count+=1
    if count==2:
        if (quadstr[0] == quadstr[3]) and (quadstr[1] == quadstr[2]):
            count = 21
    if sum([(c in abaa) for c in quadstr])!=0:
    #if '.' in quadstr:
        count = 9
    return 't%d' % count


# transfer a set of float number into ranking of [0,1]
# input a dictionary
def rank01(d):
    s = float(sum(d.values()))
    return dict((k, d[k]/s) for k in d)

# transform an np.array to [0,1] by (x-min)/(x.max - x.min)
def rankminmax(d):
	return (d-d.min())/(d.max()-d.min())    


# convert a set of float number into value of d: d = (v - mean') / std'
# a = {'a':1,'b':2,'c':3,'d':4,'e':10,'f':15}
# {'a': -0.9486, 'b': -0.6324, 'c': -0.3162, 'd': 0.0, 'e': 1.8973, 'f': 3.4785}
def rankstd(d):
    v = list(d.values())
    nv = np.array(v)
    outlier = nv.mean() + nv.std()
    #print nv.mean(), nv.std()
    bg = np.array([i for i in v if i < outlier])
    #print repr(bg)
    m = bg.mean()
    s = bg.std()
    #print m,s
    return dict((k, (d[k] - m)/s) for k in d)

# take nth root
def root(v, n):
    return v**(1./n) if 0<=v else -(-v)**(1./n)


def minpositive(d):
    return d-d.min()

# load all lines except for the empty lines
def loadlinesu(filename, ecode=None):
    import io
    with io.open(filename, encoding=ecode) as fp:
        lines = [_f for _f in (line.rstrip() for line in fp) if _f]
        return lines

# load all lines except for the empty lines
def loadlines(filename):
    with open(filename) as fp:
        lines = [_f for _f in (line.rstrip() for line in fp) if _f]
        return lines

# load lines and break each line into tuples
def loadtuples(filename, delimiter=' '):
    with open(filename) as fp:
        lines = [_f for _f in (line.rstrip() for line in fp) if _f]
        return [line.split(delimiter) for line in lines]

# load lines and break each line into tuples with regex
def loadtuplesregex(filename, regex=' '):
    import re
    with open(filename) as fp:
        lines = [_f for _f in (line.rstrip() for line in fp) if _f]
        return [re.split(regex,line) for line in lines]


# return pair tuple of (raw title and content)
def readkblockstr(linesbuff, delimiter='>'):
    raw_list = [list(x[1]) for x in groupby(linesbuff, lambda line: line[0] == delimiter)]
    return [(raw_list[i][0][1:], raw_list[i+1]) for i in range(0, len(raw_list),2)]

# convert content block into kblock format
def kblockstr(title, content, delimiter='>'):
    return '%s%s\n%s' % (delimiter, title, content)


# given two lists of coordinates. {1,..,i,..., n} in [x,y,z] format
# return RMSD
# RMSD = sqrt( 1/n * \sum_i (|| v_i - w_i ||^2)   )
def rmsd(v, w):
    if len(v) != len(w):
        print('error: vector length mismatch. v: %d w: %d' % (len(v), len(w)))
        exit(1)
    #print repr(v), repr(w)
    d = [((v[i][0]-w[i][0])*(v[i][0]-w[i][0]) + (v[i][1]-w[i][1])*(v[i][1]-w[i][1]) + (v[i][2]-w[i][2])*(v[i][2]-w[i][2])) for i in range(0, len(v))]
    #print repr(d)
    return math.sqrt(sum(d)/len(d))


# generate sm string
def smstr(npsm, alphabet):
    return '   %s\n' % ('  '.join(alphabet)) + ''.join([('%s %s\n' % (alphabet[i], ' '.join(['%2i'%n for n in npsm[i,:]]))) for i in range(len(alphabet))])

# generate sm string
def smstrf(npsm, alphabet):
	return '   %s\n' % ('  '.join(alphabet)) + ''.join([('%s %s\n' % (alphabet[i], ' '.join(['%.2f'%n for n in npsm[i,:]]))) for i in xrange(len(alphabet))])

# called in posmap_subseq()
# return a list of tuples: [(idx_long, idx_short), (, ), ...]
def subseq_align(longseq, shortseq, prefix):
    if prefix == -1:
        _err('unmatched sub sequence\nlong str: %s\nshort str: %s\n' % (longseq, shortseq))

    #print 'longseq:\n%s\n' % longseq
    #print 'shortseq:\n%s\n' % shortseq
    #print 'prefix: %d' % prefix
    # skip prefix # non-gapped alphabets
    k = 0
    while (prefix>0):
        if longseq[k] not in gaps:
            prefix-=1
        k+=1
    #print 'k: %d' % k
    retset = []
    i, j = k, 0
    while (i<len(longseq) and j<len(shortseq)):
        #print 'i:%d - %s, j:%d - %s' % (i, longseq[i], j, shortseq[j])
        if longseq[i] in gaps:
            i+=1
            continue
        if shortseq[j] in gaps:
            j+=1
            continue
        if longseq[i] == shortseq[j]:
            retset.append((i,j))
            i+=1
            j+=1
    return retset

# convert values in list into zscore
def zscore(inlist):
    n = np.array(inlist)
    return [(i-n.mean())/n.std() for i in n]


# main routine. for testing
def main():
    # test freqw
    #data = np.loadtxt('t_apply_weight.score', delimiter=',')
    data = np.loadtxt('t_lookup.score', delimiter=',')
    print(repr(data))
    w = np.loadtxt('t_apply_weight.score.70.w')
    print(repr(w))
    colset = [1,2]
    #colset = [0]
    print(repr(data[:,colset]))
    print('-----------------')
    for (c,lookup) in freqlookup(data[:,colset].T):
        print('%s%s %s %s' % (scoreaa['aa'][c[0]], scoreaa['aa'][c[1]], repr(c), ','.join([str(a) for a in lookup])))

    # test rank01
    '''
    d = {'A': 0.00005, 'B': 0.00003, 'C': 0.00001, 'D': 0.00001, 'E': 0.00000}
    d = {'A': 5, 'B': 3, 'C': 1, 'D': 1, 'E': 0}
    for k in d:
            print '%s: %.8f' % (k, d[k])
    print '-----------------'
    rd = rank01(d)
    for k in rd:
            print '%s: %.8f' % (k, rd[k])
    '''

    # test quadtype()
    '''
    qstr = ['ACAC', 'ACDG', 'ACAG', 'ACGC']
    for q in qstr:
            print (q, quadtype(q))
    '''
    # test pair_permu
    #qlist = [(('AC', 2), ('DG', 3)), (('CA', 12),('GD', 13)), (('DG',22), ('AC', 23)), (('GD', 32), ('CA', 33))]
    #qlist = [('AC','DG'), ('CA','GD'), ('DG','AC'), ('GD','CA')]
    '''
    qlist = [['A','C','D','G'], ['C','A','G','D'], ['D','G','A','C'], ['G','D','C','A']]
    for q in qlist:
            print q
            print quad_permu(q)
            print
    '''

    # test drun
    #runstr = 'ls PF*|sort' # doesn't work
    #runstr = 'pfamscan 1uhr.pdb.A.fa -json > t.json'
    #ret = drun(runstr)
    #print 'ret:\n%s\n' % ret

    # test dcall
    #ret = dcall('python utils_pfammsa.py aafreq PF00764_p90.txt')
    #print 'ret:\n%s\n' % ret
    #dcall('python utils_pfammsa.py aafreq')

    # test freq
    # print freq('jkjium')
    # defaultdict(<type 'int'>, {'i': 1, 'm': 1, 'k': 1, 'j': 2, 'u': 1})

    # test entropy
    '''
    s1 = '.........................................G.IAFSGGLDTSVAVAWMRQKG....ALPCAY.TA.........D..L...G...Q..Y.....................D......ESN...I...E.S...I.A.S....R.AK.EYG...A...E...I....A..RL.IDC..K..N.SL...VE.E.G.L..A..A............LAS..G.A...FH....IRSAGK....I....Y.FN.TTPL.G....RA.VT...G.TLLVRAML..EDN.VL.....IWG..DG.....ST.........Y.K...G..............N............D....................I...................E..R.FYRYGL.L.AN.P..........E.LKI.YKPWLDS.E.F.VAE....LGGR.KEMSD....WLKSHNLP.Y...R...D...........S........A........E........K....A..YST..DANILG.A..THE.......A....KK.LE..EL..S.....TS......IE.....I....V....E.P..............I.M......G..V....KF..W.......D...P.A....V.A.I..............T.Q..EDV..KITFKS..GR...PVAI..N......................NKD..F.SD.....PVELMKQANLIG..GRHGLG.MS.DQIENRI..IEAKS....RGI.......YEA..................PGMALLFIAYERLLSAVH.NE.E.TLA.NY.Y.Q.S.G.R.K.LGRLL.YE.........G......RW..............LDPQSL.ML.R.E.S.L.TRWV.ASA.V.SG.EVVLR.......LR...R.G.DDY.S.I.I.D.T.K.G...E...NF...S.YHPE....KLSM..E..R...T...Q.....SA.....AF.G..P..EDRIGQLTM..........RNLDIADT.....................................................................................'
    s2 = 'IAFSGGLDTSVAVAWMRQKG-ALPCAYTADLGQYDEsNIESIASRAKEYGAEIARLIDCKNSLVEEGL-AALASGAFHIrsagKIYFNTTPLGRAVTGTLLVRAMLEDNVLIWGDGSTYKGNDIERFYRYGLLANPELKIYKPWLDSEFVaelggRKEMSDWLKSHNLPYRDSAEKAYSTDANILGATHEAKKLEELSTsiEIVEPIMGVkFWDPAVAI-TQEDVKITFKSGRPVAINNKDFSDpVELMKQANLIGGRHGLGMSDQIENRIIEAKSRGIYEAPGMALLFIAYERLLSAVHNEETLANYYQSGRKLGRLLYEGRWLDPQSLMLRESLTRwVASAVSGEVVLRLRRGDdYSIIDTKGENFSYHPEKLSMERTQSaaFGPEDRIGQLT'

    s1='00113405'
    #>>> entropy("1223334444")
    #1.8464393446710154
    s1='1223334444'
    s2 = ['8.7', '9.6', '9.8', '8.8', '9.6', '10.9', '9.9', '13.9', '0.0', '10.9']
    print 'list(s1): %s' % repr(list(s1))
    b = np.array([int(i) for i in s1])
    print repr(b)
    #print 'H1: %.4f' % entropy([[float(f) for f in s2]])
    print 'H2: %.4f' % entropy([np.array([float(f) for f in s2])])
    print repr(s1)
    print 'H1: %.4f' % entropy([b])

    c = ['a','a','b','b','c']
    c = [1.2,3.4,1.2,5.5,0.0,0.0]
    print repr(c)
    print 'H(c): %.4f' % entropy([c])
    print 'H(np.array(c)): %.4f' % entropy([np.array(c)])
    '''

    #s1 = 'vLaysGGlDtsviikllkeklgeeviavavdvGqeeedldevkekalklgavksvvvDakeefvedyifpaikanalYedrYllgtalaRPliakklvevakkegaeavahGctGkGnDqvRfevsirslaPdlkviaPvRelelt....ReeeieyakekgipvevtkkkpysiDenllgrsieagiLedpknappediyeltkdpakapdepeeveiefekGvPvald....geelsv.lelieklneiagkhGvGRiDivedRlvglksReiYeapaalvLikahkdlekltlerevakfkkiveekyaelvYkGlwfsPlkealdafiektqervtGtvrvklfkGsvvvlgReseeslYdeelasydeedefdqkeaeGfikihglqakly'
    #s2 = 'vLaysGGlDtsviikllkeklgeeviavavdvGqeeedldevkekalklgavksvvvDakeefvedyifpaikanalYedrYllgtalaRPliakklvevakkegaeavahGctGkGnDqvRfevsirslaPdlkviaPvReleltReeeieyakekgipvevtkkkpysiDenllgrsieagiLedpknappediyeltkdpakapdepeeveiefekGvPvaldgeelsvlelieklneiagkhGvGRiDivedRlvglksReiYeapaalvLikahkdlekltlerevakfkkiveekyaelvYkGlwfsPlkealdafiektqervtGtvrvklfkGsvvvlgReseeslYdeelasydeedefdqkeaeGfikihglqakly'

    #s1 = 'VLAYSGGLDTSCILVWLKEQG-YDVIAYLANIGQK-EDFEEARKKALKLGAKKVFIEDVSREFVEEFIWPAIQSSALYEDRYLLGTSLARPCIARKQVEIAQREGAKYVSHGATGKGNDQVRFELSCYSLAPQIKVIAPWRMPEFYnrfkRNDLMEYAKQHGIPIPVTPKNPWSMDENLMHISYEAGILENPKNQAPPGLYTKTQDPAKAPNTPDILEIEFKKGVPVKVTnvkdGTTHQTsLELFMYLNEVAGKHGVGRIDIVENRFIGMKSRGIYETPAGTILYHAHLDIEAFTMDREVRKIKQGLGLKFAELVYTGFWHSPECEFVRHCIAKSQERVEGKVQVSVLKGQVYILGRESPLSLYNEELVSNV-QGDYEPTDATGFININSLRLKEY'
    #s2 = 'kgsvvlaysggldtscilvwlkeqgydviaylanigqkedfeearkkalklgakkvfiedvsrefveefiwpaiqssalyedryllgtslarpciarkqveiaqregakyvshgatgkgndqvrfelscyslapqikviapwrmpefynrfkrndlmeyakqhgipipvtpknpwsmdenlmhisyeagilenpknqappglytktqdpakapntpdileiefkkgvpvkvtnvkdgtthqtslelfmylnevagkhgvgridivenrfigmksrgiyetpagtilyhahldieaftmdrevrkikqglglkfaelvytgfwhspecefvrhciaksqervegkvqvsvlkgqvyilgresplslyneelvsnvqgdyeptdatgfininslrlkeyhrlqs'
    # test posmap_homoseq
    '''
    s1 = '.j..kj...'
    s2 = 'j.....k..j.'
    '''
    #retmap = posmap_homoseq(s1, s2)
    #retmap = posmap_subseq(s1,s2)
    # test posmap_subseq
    '''
    s2 = 'aa..jk.j.a'
    s1 = '.j.kj.'

    retmap = posmap_subseq(s1,s2)
    '''
    #print 's1: [%s]' % s1
    #print 's2: [%s]' % s2
    #for k,v in retmap:
    #       print 's1: %s %d -> s2: %s %d' % (s1[k], k, s2[v], v)


if __name__ == '__main__':
    main()
