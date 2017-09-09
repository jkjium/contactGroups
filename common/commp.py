import sys
import math
import operator as op

abaa = ['.','-','X','Z','U','B','O']

# AA alphabet for emboss format
aa201 = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 
		'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

aas01 = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 
		'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

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

# average resi mass, van der waals volume, polarities
aadef={
	'A':('Alanine',71.0788,67,9),
	'R':('Arginine',156.1876,148,15),
	'N':('Asparagine',114.1039,96,16),
	'D':('Aspartic.acid',115.0886,91,19),
	'C':('Cysteine',103.1448,86,7),
	'Q':('Glutamine',128.1308,114,17),
	'E':('Glutamic.acid',129.1155,109,18),
	'G':('Glycine',57.0520,48,11),
	'H':('Histidine',137.1412,118,10),
	'I':('Isoleucine',113.1595,124,1),
	'L':('Leucine',113.1595,124,3),
	'K':('Lysine',128.1742,135,20),
	'M':('Methionine',131.1986,124,5),
	'F':('Phenylalanine',147.1766,135,2),
	'P':('Proline',97.1167,90,13),
	'S':('Serine',87.0782,73,14),
	'T':('Threonine',101.1051,93,12),
	'W':('Tryptophan',186.2133,163,6),
	'Y':('Tyrosine',163.1760,141,8),
	'V':('Valine',99.1326,105,4),
	'U':('Cysteine',103.1448,86,7)
	}


# calculate n choose r
def ncr(n, r):
    r = min(r, n-r)
    if r == 0: return 1
    numer = reduce(op.mul, xrange(n, n-r, -1))
    denom = reduce(op.mul, xrange(1, r+1))
    return numer//denom


# jaccard distance for two sets
def jaccard(a, b):
	c = a.intersection(b)
	print repr(a), repr(b)
	print repr(c)
	return 1 - (float(len(c)) / (len(a) + len(b) - len(c)))

# given two lists of coordinates. {1,..,i,..., n} in [x,y,z] format
# return RMSD
# RMSD = sqrt( 1/n * \sum_i (|| v_i - w_i ||^2)   )
def rmsd(v, w):
	if len(v) != len(w):
		print 'error: vector length mismatch. v: %d w: %d' % (len(v), len(w))
		exit(1)
	#print repr(v), repr(w)
	d = [((v[i][0]-w[i][0])*(v[i][0]-w[i][0]) + (v[i][1]-w[i][1])*(v[i][1]-w[i][1]) + (v[i][2]-w[i][2])*(v[i][2]-w[i][2])) for i in xrange(0, len(v))] 
	#print repr(d)
	return math.sqrt(sum(d)/len(d))


# given two strings
# normal sequence & aligned sequence
# return map 1. key=1  pos[s1] = s2; 2. key=2 pos[s2] = s1
# s1: aligned string index, s2: pdb sequence index
def posmap(s1, s2, key=1):
	gap = ['.', '-', '_']
	ps1 = s1.translate(None, ''.join(gap))
	ps2 = s2.translate(None, ''.join(gap))
	#print 'ps1: %s\nps2: %s' % (ps1, ps2)

	retmap={}
	if ps1!=ps2:
		print 'error: not homo-str'
		print 'ps1: %s\nps2: %s' % (ps1, ps2)
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
		print 'error: incomplete map: len:%d, ps len: %d' % (len(retmap), len(ps1))
		return False

	return retmap


def getPDBUniprotMap(mapfile):
	# duplicate from utils_msa.py
	# called in utils_msa.py ncg2sdiicol()
	# load map between pdb residue ID and MSA uniprot position ID 
	# dictionary element: ('A9', (14, 'V')) : (chain+resi, (MSA position index, resn))
	# mapfile:
	# AT284 1218 T  : chain A residue T resn 284 => position 1218 (start from 0) resn T
	# AE285 1220 e  : lowercase exists!
	# AR286 -1 -	: residue number that cannot map to MSA position (does not exist)
	posmap = {}
	with open(mapfile) as fp:
		for line in fp:
			line = line.strip()
			if len(line) < 1:
				print 'getPDBUniprotMap: error loading map: %s' % mapfile
			strArray = line.split(' ')
			key = strArray[0][0] + strArray[0][2:]
			value = (int(strArray[1]), strArray[2].upper())
			posmap[key] = value
	print 'getPDBUniprotMap: %s %d maps loaded' % (mapfile, len(posmap))
	return posmap