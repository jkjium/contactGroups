import sys
import json

from utils_pfamscan import utils_pfamscan as ups
from utils_embossalign import embossalign as ea
from protein import protein

import utils_embossalign as uea
import commp as cp

class pfamscan(object):
	def __init__(self):
		pass

	# write matched HMM to as fasta file
	def dump(self):
		pass


def map_hmm2hmm(seq1, json1, seq2, json2):

	# get hmm seq
	hmm1 = json1.alnhmm.translate(None, ''.join(cp.gaps))
	hmm2 = json2.alnhmm.translate(None, ''.join(cp.gaps))
    
	# aligned position indices of two raw hmm sequences
	title = '%s-.-%s' % (json1.seqname.translate(None, ''.join(cp.illab)), json2.seqname.translate(None,''.join(cp.illab)))
	embosshmm = ea(uea.flatenalign(title, hmm1, hmm2))
	#embosshmm.dump()

	hmmmap = embosshmm.getAlignedpos()
	#print repr(hmmmap)

	# map between pfamscan hmm sequence and emboss hmm sequence
	# use emboss alignment index as key
	emboss2pfs_1 = dict((k,v) for k,v in cp.posmap_subseq(embosshmm.seqA, json1.alnhmm))
	emboss2pfs_2 = dict((k,v) for k,v in cp.posmap_subseq(embosshmm.seqB, json2.alnhmm))

	# map between pfamscan hmm 1 and pfamscan hmm 2
	pfshmm1_pfshmm2 = [(emboss2pfs_1[i], emboss2pfs_2[i]) for i in hmmmap]
	#print '-------------------------------------------------------------------------'
	#for i,j in pfshmm1_pfshmm2:
	#	print 'i:%d - %s, j:%d - %s\n' % (i, ps1.alnhmm[i], j, ps2.alnhmm[j])
	#print '-------------------------------------------------------------------------'
	#print '-------------------------------------------------------------------------'
	#for i,j in pfshmm1_pfshmm2:
	#	print 'i:%d - %s, j:%d - %s\n' % (i, ps1.alnseq[i], j, ps2.alnseq[j])
	#print '-------------------------------------------------------------------------'

	#print 'mapped pfshmm1: %s\n' % ''.join([ps1.alnhmm[p[0]] for p in pfshmm1_pfshmm2])
	#print 'mapped pfshmm2: %s\n' % ''.join([ps2.alnhmm[p[1]] for p in pfshmm1_pfshmm2])

	seq2alnhmm_1 = dict((k,v) for k,v in cp.posmap_subseq(json1.alnseq, seq1))
	seq2alnhmm_2 = dict((k,v) for k,v in cp.posmap_subseq(json2.alnseq, seq2))
	#print '-------------------------------------------------------------------------'
	#for k in seq2alnhmm_1:
	#	print 'k:%d - %s, v:%d - %s\n' % (k, ps1.alnseq[k], seq2alnhmm_1[k], s1[seq2alnhmm_1[k]])
	#print '-------------------------------------------------------------------------'
	'''
	print repr(pfshmm1_pfshmm2)
	print 'ps1.alnhmm:\n%s\n' % ps1.alnhmm
	print 'ps1.alnseq:\n%s\n' % ps1.alnseq
	print 'ps2.alnhmm:\n%s\n' % ps2.alnhmm
	print 'ps2.alnseq:\n%s\n' % ps2.alnseq
	'''
	# hmm index may not in alnseq index because alnseq are gapped in order to align to hmm profile
	seq_map = [(seq2alnhmm_1[i], seq2alnhmm_2[j]) for i,j in pfshmm1_pfshmm2 if (i in seq2alnhmm_1 and j in seq2alnhmm_2)]
	'''
	print 'mapped seq1:\n%s' % ''.join(s1[p[0]] for p in seq_map)
	print 'mapped seq2:\n%s' % ''.join(s2[p[1]] for p in seq_map)
	'''
	'''
	for k,v in seq_map:
		print 's1: %s %d -> s2: %s %d' % (seq1[k], k, seq2[v], v)
	'''
	return seq_map


# map pdb resi to msa index
def pdbResi2MSA(pdbfile, pdbseqfafile, pdbjsonfile, msafafile, msajsonfile, pfamid):
	p = protein(pdbfile)

	# load sequences
	pdbseq = [s for s in cp.fasta_iter(pdbseqfafile)][0][1]
	msaseq = [s for s in cp.fasta_iter(msafafile)][0][1]

	# load pfamscan json object
	pdbjson = ups(pdbjsonfile).getMatchpfs(pfamid)
	msajson = ups(msajsonfile).getMatchpfs(pfamid)

	# get map between pdb pos and msa pos
	pdbpos2msapos = dict((k,v) for k,v in map_hmm2hmm(pdbseq, pdbjson, msaseq, msajson))

	# replace pdb pos with pdb resi
	resi2msa = [(p.ca[i], pdbpos2msapos[i]) for i in xrange(0, len(p.ca)) if i in pdbpos2msapos]

	outfile = '%s-%s.map' % (pdbfile, pfamid)
	outstr = '\n'.join(['%d %s %d %s' % (k.resSeq, cp.aa2a[k.resName], v, msaseq[v]) for (k,v) in resi2msa])
	with open(outfile, 'w') as fp:
		fp.write(outstr)
	print 'save to %s' % outfile
	'''
	for k,v in resi2msa:
		print 'pdb %d - %s, msa %d - %s' % (k.resSeq, cp.aa2a[k.resName], v, msaseq[v])
	'''

	return outstr	

def resi2msa():
	if len(sys.argv) < 8:
		print 'Usage: python utils_resimap.py resi2msa pdbfile pdbseqfafile pdbjsonfile msafafile msajsonfile pfamid'
		print '$ python utils_resimap.py resi2msa 1ni3.pdb PF06071_pdb.fa PF06071_1ni3.json PF06071_MSA.fa PF06071.json PF06071'
		exit()
	print pdbResi2MSA(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7])


def test():
	# 1ni3.pdb: 			raw pdb
	# PF06071_pdb.fa: 		pdb seq	
	# PF06071_1ni3.json: 	pfamscan result from PF06071_pdb.fa
	# PF06071_MSA.fa: 		raw MSA sequence (with gap)
	# PF06071.json: 		pfamscan result from ungapped PF06071_seq.fa
	pdbResi2MSA('1ni3.pdb', 'PF06071_pdb.fa', 'PF06071_1ni3.json', 'PF06071_MSA.fa', 'PF06071.json', 'PF06071')
	pass

# main routine
def main():
	if len(sys.argv)<2:
		print 'Usage: python utils_protein.py cmd pdbfile [args ...]'
		return

	dispatch = {
		'test':test,
		'resi2msa':resi2msa
	}

	cmd = sys.argv[1]

	if sys.argv[1] not in dispatch:
		print 'invalid cmd: %s' % sys.argv[1]
		return
	else:
		dispatch[sys.argv[1]]()

if __name__ == '__main__':
	main()