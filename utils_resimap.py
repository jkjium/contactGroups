import sys
import json
import collections

from utils_pfamscan import utils_pfamscan as ups
from utils_embossalign import embossalign as ea
from protein import protein
from utils_testcase import palign

import utils_embossalign as uea
import commp as cp
import numpy as np

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

	match_level, hmmmap = embosshmm.getAlignedpos()
	cp._info('info:%s hmm match level: %.2f' % (title, match_level))
	if match_level < 0.9:
		cp._info('err:%s hmm match less than 90% : %.2f' % (title, match_level))
		exit()

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
# pdbfile: 		pdb structure file
# pdbseqfafile: get from utils_protein writeseqfa pdbfile {chain}
# pdbjsonfile: 	get from pfamscan pdbseqfafile
# msafafile:	MSA sequence WITH GAPs extracted from pfam MSA
# msajsonfile:	get from pfamscan MSA sequence WITHOUT GAPs
# output: resi resn msai msan
def pdbResi2MSA(pdbfile, chainid, pdbseqfafile, pdbjsonfile, msafafile, msajsonfile, pfamid):
	p = protein(pdbfile, chain=chainid)

	# load sequences
	pdbseq = [s for s in cp.fasta_iter(pdbseqfafile)][0][1]
	msa1 = [s for s in cp.fasta_iter(msafafile)][0]
	msaseq = msa1[1]
	fullhead = msa1[0]
	strarr = fullhead.split('/')
	msahead = strarr[0]

	# load pfamscan json object
	pdbjson = ups(pdbjsonfile).getMatchpfs(pfamid)
	if pdbjson == False:
		cp._info('err: %s not found in %s' % (pfamid, pdbjsonfile))
		return
	msajson = ups(msajsonfile).getMatchpfs(pfamid)
	if msajson == False:
		cp._info('err: %s not found in %s' % (pfamid, msajsonfile))
		return

	# get map between pdb pos and msa pos
	pdbpos2msapos = dict((k,v) for k,v in map_hmm2hmm(pdbseq, pdbjson, msaseq, msajson))

	# replace pdb pos with pdb resi
	#resi2msa = [(p.ca[i], pdbpos2msapos[i]) for i in xrange(0, len(p.ca)) if i in pdbpos2msapos]
	# resArray, gives a list of keys eg. (A,Q,70), (A,I,71), (A,V,72)
	resi2msa = [(p.resArray[i], pdbpos2msapos[i]) for i in xrange(0, len(p.resArray)) if i in pdbpos2msapos]

	outfile = '%s-%s-%s.map' % (pdbfile, chainid, pfamid)
	#outstr = '\n'.join(['%d %s %d %s' % (k.resSeq, cp.aa2a[k.resName], v, msaseq[v]) for (k,v) in resi2msa])
	outstr = '\n'.join(['%d %s %d %s' % (k[2], k[1], v, msaseq[v]) for (k,v) in resi2msa])
	with open(outfile, 'w') as fp:
		fp.write(outstr)
	cp._info('save map to %s' % outfile)
	'''
	for k,v in resi2msa:
		print 'pdb %d - %s, msa %d - %s' % (k.resSeq, cp.aa2a[k.resName], v, msaseq[v])
	'''

	return outstr	

# main mapping function
def resi2msa(arglist):
	if len(arglist) < 7:
		print 'Usage: python utils_resimap.py resi2msa pdbfile pdbseqfafile pdbjsonfile msafafile msajsonfile pfamid'
		print '$ python utils_resimap.py resi2msa 1a8v.pdb B 1a8v.pdb.B.fa 1a8v.pdb.B.fa.json PF07497_p90_MSA.fa PF07497_p90_seq.fa.json PF07497'
		exit()
	#pdbResi2MSA(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7])
	pdbResi2MSA(arglist[0], arglist[1], arglist[2], arglist[3], arglist[4], arglist[5], arglist[6])


# given a sequence generate Pfam MSA alignment
# output: .fa file that aligns to the existing MSA
# seq.fa
# seq.fa.json : pfamscan seq.fa > seq.fa.json
# PF00000_MSA(seq).fa : python utils_embossalign.py findsimilar seq.fa PF00000_full.txt  > seq-PF00000.report &
#                       python utils_pfammsa.py getsinglemsa PF00000_full.txt A0A1S3GL90_DIPOR/587-617 PF00000
# PF00000_seq.fa.json: pfamscan PF00000_seq.fa > PF00000_seq.fa.json
def appendseq2msa(arglist):
	if len(arglist) < 5:
		cp._info('Usage: python utils_resimap.py appendseq2msa seq.fa seq.fa.json PF00000_MSA.fa PF00000_seq.fa.json pfamid')

	seqfafile = arglist[0]
	seqjsonfile = arglist[1]
	msafafile = arglist[2]
	msajsonfile = arglist[3]
	pfamid = arglist[4]

	# load sequences
	seq = cp.fasta_iter(seqfafile).next()[1]
	msahead, msaseq = cp.fasta_iter(msafafile).next()

	# load pfamscan json object
	seqjson = ups(seqjsonfile).getMatchpfs(pfamid)
	if seqjson == False:
		cp._err('%s not found in %s' % (pfamid, seqjsonfile))
	msajson = ups(msajsonfile).getMatchpfs(pfamid)
	if msajson == False:
		cp._err('%s not found in %s' % (pfamid, msajsonfile))

	# get map between pdb pos and msa pos
	# v is msapos, k is seqpos
	# different from resi2msa
	msapos2seqpos = dict((v,k) for k,v in map_hmm2hmm(seq, seqjson, msaseq, msajson))
	outseq = [] 
	for i in xrange(0, len(msaseq)):
		if i in msapos2seqpos:
			outseq.append(seq[msapos2seqpos[i]])
		else:
			outseq.append('.')

	outfa = '>%s\n%s\n' % (seqfafile, ''.join(outseq))
	outfile = '%s-%s.fa' % (seqfafile, pfamid)
	with open(outfile, 'w') as fout:
		fout.write(outfa.upper())
	cp._info('save msa to %s' % (outfile))

	outfile = '%s-%s.map' % (seqfafile, pfamid)
	klist = msapos2seqpos.keys()
	klist.sort()
	# format: resi resn msai msan
	outstr = '\n'.join(['%d %s %d %s' % (msapos2seqpos[k], seq[msapos2seqpos[k]], k, msaseq[k]) for k in klist])
	with open(outfile, 'w') as fout:
		fout.write(outstr.upper())
	cp._info('save map to %s' % outfile)


# generate an (aligned)MSA sequence from the resimapfile
# input: resimapfile, output_seq_header, template_MSA_seq(FASTA), outputfile
# output: the generated MSA sequence in FASTA format
def createmsaseq(arglist):
	if len(arglist) < 4:
		cp._err('Usage: python utils_resimap.py createmsaseq 6v5d_0.pdb-A-PF00240.map header PF00240_MSA.fa outfile')

	mapfile = arglist[0]
	outheader = arglist[1]
	msafile = arglist[2]
	outfile = arglist[3]

	# load resimap
	# ri rn msai msan
	# 3 V 39 I
	# msai2resn[39] = 'V'
	msai2resn = {}
	for line in cp.loadlines(mapfile):
		sarr = line.split(' ')
		msai2resn[int(sarr[2])] = sarr[1]

	#print len(msai2resn)
	#print msai2resn

	# read template MSA seq
	templatehead, templateseq = cp.fasta_iter(msafile).next()

	outseq = [] 
	for i in xrange(0, len(templateseq)):
		if i in msai2resn:
			outseq.append(msai2resn[i])
		else:
			outseq.append('.')

	outfa = '>%s\n%s\n' % (outheader, ''.join(outseq))
	with open(outfile, 'w') as fout:
		fout.write(outfa)

	#print '>%s\n%s\n' % (templatehead, templateseq)
	#print outfa
	cp._info('save to %s' % (outfile))



# append resi to tsvfile according the column id in tsvfile
def appendresi(arglist):
	if len(arglist) < 4:
		cp._err('Usage: python utils_resimap.py appendresi tsvfile column_index mapfile outfile')

	infile = arglist[0]
	col = int(arglist[1])
	mapfile = arglist[2]
	outfile = arglist[3]

	residict = collections.defaultdict(lambda:-1)
	with open(mapfile) as fp:
		for line in fp:
			line = line.strip()
			if len(line)==0:
				continue
			# resi resn msai msan
			# 406 K 269 G
			sarr = line.split(' ')
			residict[sarr[2]] = sarr[0]
	#cp._info('%d map loaded.' % len(residict))

	fout = open(outfile, 'w')
	# load input file
	with open(infile) as fp:
		for line in fp:
			line = line.strip()
			if len(line) == 0:
				continue
			sarr = line.split(' ')
			fout.write('%s %s\n' % (line, residict[sarr[col]]))
	fout.close()
	cp._info('save to %s' % outfile)	 	


def cg2msai_sdii(arglist):
	if len(arglist) < 5:
		cp._err('Usage: python utils_resimap cg2msai cgfile mapfile_left mapfile_right sdiifile outfile')

	cgfile = arglist[0]
	mapfile1 = arglist[1]
	mapfile2 = arglist[2]
	sdiifile = arglist[3]
	outfile = arglist[4]

	# load map
	resmap1 = {}
	with open(mapfile1) as fp:
		for line in fp:
			line = line.strip()
			if len(line) == 0:
				continue
			sarr = line.split(' ')
			#resi resn msai msan
			#149 P 88 Q
			#resmap[sarr[2]] = sarr[0]
			#		resi 		msai
			resmap1[sarr[0]] = sarr[2]


	resmap2 = {}
	with open(mapfile2) as fp:
		for line in fp:
			line = line.strip()
			if len(line) == 0:
				continue
			sarr = line.split(' ')
			#resi resn msai msan
			#149 P 88 Q
			#resmap[sarr[2]] = sarr[0]
			#		resi 		msai
			resmap2[sarr[0]] = sarr[2]

	# load sdii
	sdii = {}
	with open(sdiifile) as fp:
		for line in fp:
			line = line.strip()
			if len(line) == 0:
				continue
			sarr = line.split(' ')
			#msai msai mi resi1 resi1 resi2 resi2 .... 
			#326 653 0.315745486152245 232 328 117 213 115 211
			sdii['%s %s' % (sarr[0], sarr[1])] = sarr[2]

	#print repr(sdii)
	outstr = []
	with open(cgfile) as fp:
		for line in fp:
			line = line.strip()
			if len(line) == 0:
				continue
			sarr = line.split(' ')
			#  map1 map2
			#122 F 254 W chain E A
			if (sarr[0] in resmap1) and (sarr[2] in resmap2):
				msai1 = resmap1[sarr[0]]
				msai2 = resmap2[sarr[2]]
				sdiikey = '%s %s' % (msai2, msai1) if int(msai1) >= int(msai2) else '%s %s' % (msai1, msai2)
				#print sdiikey
				if sdiikey in sdii:
					outstr.append('%s %s %s %s' %(line, msai1, msai2, sdii[sdiikey]))

	with open(outfile, 'w') as fout:
		fout.write('%s\n' % '\n'.join(outstr))
	cp._info('save %d records to %s' % (len(outstr), outfile))	

# read index from .vec (single column) file
# map loaded indices according to flag: resi2msai
# save to a .vec file 
def mapvec(arglist):
	if len(arglist) < 4:
		cp._err('Usage: python utils_resimap.py mapvec resi2msai|msai2resi mapfile col.vec outfile')
	flag = arglist[0]
	mapfile = arglist[1]
	vecfile = arglist[2]
	outfile = arglist[3]

	# 31 K 56 R
	mlist = []
	for line in cp.loadlines(mapfile):
		sarr = line.split(' ')
		mlist.append((sarr[0], sarr[2]))

	mdict = {}
	if flag == 'resi2msai':
		# given resi output mesi
		for m in mlist:
			mdict[m[0]] = m[1]
	elif flag == 'msai2resi':
		# given msai output resi
		for m in mlist:
			mdict[m[1]] = m[0]
	else:
		cp._err('Wrong flag : %s' % flag)
	
	outlist=[]
	for line in cp.loadlines(vecfile):
		if line in mdict:
			outlist.append(mdict[line])
		else:
			outlist.append('-191')
	with open(outfile, 'w') as fout:
		fout.write('%s\n' % ('\n'.join(outlist)))

	cp._info('save to %s' % outfile)



# column (group) to resi
# convert column (group) into resi
# one group per line
def col2resi(arglist):
	if len(arglist) < 3:
		cp._err('Usage: python utils_resimap.py col2resi mapfile colfile outresifile')

	mapfile = arglist[0]
	colfile = arglist[1]
	outfile = arglist[2]

	resmap = {}
	with open(mapfile) as fp:
		for line in fp:
			line = line.strip()
			if len(line) == 0:
				continue
			sarr = line.split(' ')
			#resi resn msai msan
			#149 P 88 Q
			resmap[sarr[2]] = sarr[0]

	outstr = []
	with open(colfile) as fp:
		for line in fp:
			line = line.strip()
			if len(line) == 0:
				continue
			sarr = line.split(' ')
			#326 653 
			outstr.append(' '.join([resmap[c] for c in sarr]))

	with open(outfile, 'w') as fout:
		fout.write("%s\n" % ('\n'.join(outstr)))
	cp._info('save to %s' % outfile)

# reverse of col2resi
# map a group of resi to msai
def resi2col(arglist):
	if len(arglist) < 3:
		cp._err('Usage: python utils_resimap.py resi2col mapfile resifile outfile')
	mapfile = arglist[0]
	resifile = arglist[1]
	outfile = arglist[2]

	resmap={}
	# resi resn msai msan
	# 101 K 2 K
	for line in cp.loadlines(mapfile):
		sarr = line.split(' ')
		resmap[sarr[0]] = sarr[2]

	outstr = []
	for line in cp.loadlines(resifile):
		sarr = line.split(' ')
		# neighbors from pdb file may not be included in the map file
		cols = []
		for r in sarr:
			if r not in resmap:
				#continue
				cols.append('-1')
			else:
				cols.append(resmap[r])
		outstr.append(' '.join(cols))

	with open(outfile, 'w') as fout:
		fout.write("%s\n" % ('\n'.join(outstr)))
	cp._info('save to %s' % outfile)

# append msai to dca output
# for cecolumn, cflat
# $ python utils_resimap.py dca2msa PF00098_full.txt PF00098_full.txt.dca PF00098_full.rdca
def dca2msa(arglist):
	if len(arglist) < 4:
		cp._err('Usage: python utils_resimap.py dca2msa msafile seqheader dcafile outfile')

	msafile = arglist[0]
	msahead = arglist[1]
	dcafile = arglist[2]
	outfile = arglist[3]

	# get the specified (msahead) entry
	seq = ''
	msa = ''
	resi_start = -1
	for h, s in cp.fasta_iter(msafile):
		if msahead == h:
			msa = s
			seq = s.translate(None, ''.join(cp.gaps))
			resi=h.split('/')[1]
			resi_start = int(resi.split('-')[0])

	if resi_start == -1:
		cp._err('%s not found' % msahead)

	'''
	head, msa = cp.fasta_iter(msafile).next()
	seq = msa.translate(None, ''.join(cp.gaps))
	#print '%s\n%s\n%s' % (head, msa, seq)
	resi=head.split('/')[1]
	resi_start = int(resi.split('-')[0])
	'''

	# seq2msa = dict((k+resi_start,v) for k,v in cp.posmap_subseq(seq, msa))
	seq2msa = dict((k,v) for k,v in cp.posmap_subseq(seq, msa))
	# dca index = k+resi_start 
	#for k in seq2msa:
	#	print '%d %d %s %d' % (k, seq2msa[k], msa[seq2msa[k]], k+resi_start)

	fout = open(outfile ,'w')
	# convert dca index
	with open(dcafile) as fp:
		for line in fp:
			line = line.strip()
			if len(line) == 0:
				continue
			# 240 V 241 E 0.424638 0.170393
			sarr = line.split(' ')
			fout.write('%d %d %s\n' % (seq2msa[int(sarr[0])-resi_start], seq2msa[int(sarr[2])-resi_start], line))
	fout.close()
	cp._info('save to %s' % outfile)

# after dca2msa => .rdca
# for the hybrid procedure
# read .rcol from msareduce
# map column id to .rdca file
def dcamsa2rcol(arglist):
	if len(arglist) < 3:
		cp._err('Usage: python utils_resimap.py dca2msabyrcol .rdcafile .rcolfile .cdcaoutfile')

	rdcafile = arglist[0]
	rcollist = [int(i) for i in np.loadtxt(arglist[1], delimiter=',')]
	outfile = arglist[2]
	
	fout = open(outfile, 'w')
	# line: 0 51 0 N 5 P 0.073385 0.000022
	for line in cp.loadlines(rdcafile):
		sarr = line.split(' ')

		idx1 = int(sarr[0])
		idx2 = int(sarr[1])

		col1 = rcollist[idx1]
		col2 = rcollist[idx2]

		fout.write('%d %d %s\n' % (col1, col2, line))

	fout.close()
	cp._info('save to %s' % outfile)
	


# msa2pdb: resolve pdb seq and Pfam MSA pdb seq mapping problem (for the same protein described in pdb structure)
# 1. full pdb sequence is longer than pfam domain
# 2. domain sequence in pdb and in the MSA are different
# use local alignment to get correct alignment
# for p53/rmsd_comparison
# input: msa.seq(gaps are allowed), pdb.seq, pdb_msa.water.align
# output: msa_pdb.map (given msai output pdbi)
def msa2pdb(arglist):
	if len(arglist) < 4:
		cp._err('Usage:python utils_resimap.py msa2pdb msa.seq pdb.seq msa_pdb.align.flat outfile')
	msaseqfile = arglist[0]
	pdbseqfile = arglist[1]
	flatfile = arglist[2]
	outfile = arglist[3]

	msaseq = cp.loadlines(msaseqfile)[0]
	#print 'msa sequence:'
	#print msaseq
	pdbseq = cp.loadlines(pdbseqfile)[0]
	#print 'pdb sequence:'
	#print pdbseq
	flatstr = cp.loadlines(flatfile)[0]
	#print 'flat content:'
	#print flatstr
	pa =  palign(flatstr)
	#pa.dump()
	alnstr_msa = pa.seqA
	#print 'msa aln seq:'
	#print alnstr_msa
	alnstr_pdb = pa.seqB
	#print 'pdb aln seq:'
	#print alnstr_pdb
	apos = pa.alnpos()
	#print apos

	# given aligned msa index, output msai
	map_msa2aln = cp.posmap_subseq(alnstr_msa, msaseq)
	#print map_msa2aln
	#print msaseq[520], alnstr_msa[177]
	# given aligned pdb index, output pdbi
	map_aln2pdb = cp.posmap_subseq(alnstr_pdb, pdbseq)
	#print map_aln2pdb
	#print alnstr_pdb[163], pdbseq[157]

	# combine result into a map, msai to pdbi
	#map_msa2pdb = dict((map_msa2aln[i], map_aln2pdb[i]) for i in apos)
	map_msa2pdb = [(map_msa2aln[i], map_aln2pdb[i]) for i in apos]
	#print map_msa2pdb
	#print msaseq[224], pdbseq[]

	with open(outfile, 'w') as fout:
		fout.write('\n'.join(['%d %d' % (i[0], i[1]) for i in map_msa2pdb]))
	cp._info('save msa2pdb map to %s' % outfile)

# given a column of msai output the resi
# input infile, msai_index, mapfile, outfile
# print mapped resi
def msa2rescol(arglist):
	if len(arglist) < 3:
		cp._err('Usage: python utils_resimap.py msa2rescol PF00418.rdca 0,1 2mz7_0.pdb-A-PF00418.map')
	infile = arglist[0]
	idx = [int(i) for i in arglist[1].split(',')]
	mapfile = arglist[2]

	# 2mz7_0.pdb-A-PF00418.map
	# res n msa n
	# 275 V 14 V
	resimap = {}
	for line in cp.loadlines(mapfile):
		sarr= line.split(' ')
		# given msai output resi
		resimap[sarr[2]] = sarr[0]

	# load msai using infile and idx
	#==> PF00418.rdca
	#14 15 587 V 588 Q 0.30887 0.2327
	outresi = []
	for line in cp.loadlines(infile):
		sarr = line.split(' ')
		reslist = []
		for i in idx:
			if sarr[i] in resimap:
				reslist.append('%s' % resimap[sarr[i]])
			else:
				reslist.append('-191')
		print ' '.join(reslist)





# append resi in front of the result from mp_ce_sdii (weight)
def sdii2resi(arglist):
	if len(arglist) < 3:
		cp._err('Usage: python utils_resimap sdii2res sdiifile mapfile outfile')

	sdiifile = arglist[0]
	mapfile = arglist[1]
	outfile = arglist[2]

	# load map
	resmap = {}
	with open(mapfile) as fp:
		for line in fp:
			line = line.strip()
			if len(line) == 0:
				continue
			sarr = line.split(' ')
			#resi resn msai msan
			#149 P 88 Q
			resmap[sarr[2]] = sarr[0]

	outstr = []
	with open(sdiifile) as fp:
		for line in fp:
			line = line.strip()
			if len(line) == 0:
				continue
			sarr = line.split(' ')
			#326 653 0.315745486152245
			if (sarr[0] in resmap) and (sarr[1] in resmap):
				outstr.append('%s %s %s' %(line, resmap[sarr[0]], resmap[sarr[1]]))

	with open(outfile, 'w') as fout:
		fout.write('%s\n' % '\n'.join(outstr))
	cp._info('save %d records to %s' % (len(outstr), outfile))



# append resi to triplets SDII output 
def triplet2resi(arglist):
	if len(arglist) < 3:
		cp._err('Usage: python utils_resimap triplet2resi sdiifile mapfile outfile')

	sdiifile = arglist[0]
	mapfile = arglist[1]
	outfile = arglist[2]

	# load map
	resmap = {}
	with open(mapfile) as fp:
		for line in fp:
			line = line.strip()
			if len(line) == 0:
				continue
			sarr = line.split(' ')
			#resi resn msai msan
			#149 P 88 Q
			resmap[sarr[2]] = sarr[0]

	outstr = []
	with open(sdiifile) as fp:
		for line in fp:
			line = line.strip()
			if len(line) == 0:
				continue
			sarr = line.split(' ')
			#326 653 0.315745486152245
			if (sarr[0] in resmap) and (sarr[1] in resmap) and (sarr[2] in resmap):
				outstr.append('%s %s %s %s' % (line, resmap[sarr[0]], resmap[sarr[1]], resmap[sarr[2]]))

	with open(outfile, 'w') as fout:
		fout.write('%s\n' % '\n'.join(outstr))
	cp._info('save %d records to %s' % (len(outstr), outfile))









'''
def test():
	# 1ni3.pdb: 			raw pdb
	# PF06071_pdb.fa: 		pdb seq	
	# PF06071_1ni3.json: 	pfamscan result from PF06071_pdb.fa
	# PF06071_MSA.fa: 		raw MSA sequence (with gap)
	# PF06071.json: 		pfamscan result from ungapped PF06071_seq.fa
	pdbResi2MSA('1ni3.pdb', 'PF06071_pdb.fa', 'PF06071_1ni3.json', 'PF06071_MSA.fa', 'PF06071.json', 'PF06071')
	pass

'''
# main routine
'''
def main():
	if len(sys.argv)<2:
		print 'Usage: python utils_protein.py cmd pdbfile [args ...]'
		return

	dispatch = {
		'test':test,
		'resi2msa':resi2msa,
		'dca2msa': dca2msa
	}

	if sys.argv[1] not in dispatch:
		print 'invalid cmd: %s' % sys.argv[1]
		return
	else:
		dispatch[sys.argv[1]](sys.argv[2:])
'''

if __name__ == '__main__':
	cp.dispatch(__name__)
