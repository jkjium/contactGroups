import sys
import numpy as np

from atom import atom
from msa import msa
from protein import protein

import common.commp as cp

def addscore(msacol, sm):
	for aa in msacol:
		#if aa.upper() not in cp.abaa:
		if aa.upper() in cp.aa201:
			sm[aa.upper()]+=1

# check whether there is unmatched residue in the map file
def validncg(ncg, rtmap):
	for a in ncg:
		if a.resName.upper() in cp.abaa or a.resSeq <0:
			return False
		k = '%s%d' % (a.chainID, a.resSeq)
		if rtmap[k][0] == -1:
			return False
	return True

# a chuck of output produced from a pfam msa 
def cgmsa_composition(pdbfile, msafile, mapfile):
	'''
	find corresponding columns in msa
	'''
	cgsize = 2 # pairwised contact 
	# init score

	m = msa(msafile)
	# matrix format of msa
	# for columnwise access
	msamat = np.array([list(s[1]) for s in m.msaArray]) 

	# 'A9': (14, 'V')
	rtmap = cp.getPDBUniprotMap(mapfile)

	p = protein(pdbfile)

	ret = []
	# each ncg yields a result
	for ncg in p.contactbynearest(p.atoms,cgsize):
		if validncg(ncg, rtmap) == False:
			#print 'skip invalid ncg: %s' % ' '.join(['%s%d' % (a.chainID, a.resSeq) for a in ncg])
			continue

		pdbcg = ['%s%d' % (a.chainID, a.resSeq) for a in ncg]
		pdbcgstr = ','.join(pdbcg)

		mass = sum([cp.aadef[cp.aa2a[a.resName]][1] for a in ncg])
		volume = sum([cp.aadef[cp.aa2a[a.resName]][2] for a in ncg])
		polar = sum([cp.aadef[cp.aa2a[a.resName]][3] for a in ncg])
		# a list (two for pairwise) of residues in the contact. e.g. (A6,A112)

		msacg = [rtmap[k][0] for k in pdbcg]
		msacgstr = ','.join([str(col) for col in msacg])

		sm = {}
		for aa in cp.aa201:
			sm[aa]=0

		for col in msacg:
			addscore(list(msamat[:,col]), sm)

			#print sm[aa], 2*len(m.msaArray)
		compstr = ','.join([str(int(5.0*sm[aa]/len(m.msaArray))) for aa in cp.aa201])

		ret.append('%s %s %s %s %s %.2f %.2f %.2f' % (pdbfile, msafile, pdbcgstr, msacgstr, compstr, mass, volume, polar))
		#print '%s %s %s %s %s %.2f %.2f %.2f\n' % (pdbfile, msafile, pdbcgstr, msacgstr, compstr, mass, volume, polar)

	return '\n'.join(ret)


def main():
	'''
	For residue contact msa column composition
	covariate data 
	read pdb, msa
	output pdb+pfam+composition
	'''

	if len(sys.argv) < 2:
		print 'Usage: python proc_cv_composition 13-talign-ok.tsv'
		return

	blackboard = sys.argv[1]
	outfile = 'covariate.%s.rpcomp' % blackboard
	fout = open(outfile, 'w')
	with open(blackboard) as fp:
		for line in fp:
			strArray = line.strip().split(' ')
			pdbfile = '%s-%s.rpdb.tip' % (strArray[0].lower(), strArray[1])
			pfamfile = '%s_full.txt.rseq' % strArray[3] 
			rtmapfile = '%s-%s-%s-%s.map' % (strArray[0].lower(), strArray[1], strArray[3], strArray[6])
			#1a0p-A-PF00589-XERD_ECOLI.map
			fout.write(cgmsa_composition(pdbfile, pfamfile, rtmapfile))
	fout.close()
	print 'save to %s.' % outfile

if __name__ == '__main__':
	main()