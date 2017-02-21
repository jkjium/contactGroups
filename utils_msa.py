'''
get msa position id 
'''
import sys
import numpy as np
import collections
import math
from msa import msa
from protein import protein

# print out whether a uniprot name is existing in the PFam MSA
def printUniprot():
	if len(sys.argv) < 4:
		print 'showUniprot: print existing target uniprot in PFXXXXX'
		print 'example: python utils_msa.py printuniprot PF00154_full.txt Q7X416'
		print 'output: PF00154_full.txt Q7X416 1'
		return

	msafile = sys.argv[2]
	target = sys.argv[3]

	m = msa(msafile)
	outArray = m.searchUniprot(target)
	if len(outArray) == 0:
		print '%s %s 0' % (msafile, target)
	else:
		for name in outArray:
			print '%s %s %d' % (msafile, name, len(outArray), ' '.join(outArray))



def sdii2resi():
	if len(sys.argv) < 5:
		print 'resi2target: given a residue number output the corresponding position in target msa'
		print 'example:python utils_msa.py sdii2resi PF07714_full.fa.r50 BTK_HUMAN 1k2p.pdb PF07714_full.fa.r50.3128_3_sdii\n'
		print 'output: PF07714_full.fa.r50.3128_3_sdii_resi'
		return

	msafile = sys.argv[2]
	target = sys.argv[3]
	pdbfile = sys.argv[4]
	sdiifile = sys.argv[5]

	print 'msafile: %s\ntarget header: %s\npdbfile: %s\nsdii file: %s' % (msafile, target, pdbfile, sdiifile)
	m = msa(msafile)
	p = protein(pdbfile)
	rtmap = m.getResiTargetMap(p, target)
	if len(rtmap) < 1:
		print 'error occoured in generating rtmap'
		return
	#print '%s: %s' % (tvar, repr(rtmap[tvar]))
	# construct trmap from rtmap
	# 3128: (B641, 'R')
	trmap = {}
	#trmap = {v: k for k, v in rtmap.iteritems()}
	for k in rtmap:
		msai, resn = rtmap[k]
		if msai in trmap:
			print 'error. duplicate key [%d] in rtmap' % msai
			return
		trmap[msai] = (k, resn)

	#print trmap

	# read sdii file
	with open(sdiifile) as f:
		sdiilines = f.readlines()

	outfile = sdiifile + '_resi'
	fout = open(outfile, 'w')

	# 52 [pid:20029] 926-3089-3128 0.001106226720675
	count = 0
	for line in sdiilines:
		count += 1
		print '%d/%d processed ...' % (count, len(sdiilines))
		strArr = line.strip().split(' ')
		msailist = strArr[2].split('-')
		sdiivalue = strArr[3]
		fout.write('%s %s\n' % ('-'.join([repr(trmap[int(i)]) for i in msailist]), sdiivalue))
	fout.close()
	print 'done.\noutput file: [%s]' % outfile


def msai2resi():
	if len(sys.argv) < 4:
		print 'msai2resi: output the mapping between msa position index and pdb residue number'
		print 'example:python utils_msa.py msai2resi PF07714_full.fa BTK_HUMAN 1k2p.pdb\n'
		print 'output: PF07714_full.fa.1k2p.pdb.map'
		return

	msafile = sys.argv[2]
	target = sys.argv[3]
	pdbfile = sys.argv[4]
	outfile = msafile+'.'+pdbfile+'.map'

	print 'msafile: %s\ntarget header: %s\npdbfile: %s\noutput file: %s' % (msafile, target, pdbfile, outfile)
	m = msa(msafile)
	p = protein(pdbfile)
	rtmap = m.getResiTargetMap(p, target)
	if len(rtmap) < 1:
		print 'error occoured in generating rtmap'
		return
	#print '%s: %s' % (tvar, repr(rtmap[tvar]))
	# construct trmap from rtmap
	# 3128: (B641, 'R')
	trmap = {}
	#trmap = {v: k for k, v in rtmap.iteritems()}
	fout = open(outfile ,'w')
	for k in rtmap:
		msai, resn = rtmap[k]
		if msai in trmap:
			print 'error. duplicate key [%d] in rtmap' % msai
			return
		trmap[msai] = (k, resn)
		fout.write('%d %d %d' % (msai, k, resn))
	fout.close()
	#print trmap	

# output: B641: (3128, 'R')
def resi2target():
	if len(sys.argv) < 5:
		print 'resi2target: given a residue number output the corresponding position in target msa'
		print 'example:python utils_msa.py resi2target PF07714_full.fa.r50 BTK_HUMAN 1k2p.pdb B641\n'
		return

	msafile = sys.argv[2]
	target = sys.argv[3]
	pdbfile = sys.argv[4]
	tvar = sys.argv[5]

	print 'msafile: %s\ntarget header: %s\npdbfile: %s\ntarget variable: %s' % (msafile, target, pdbfile, tvar)
	m = msa(msafile)
	p = protein(pdbfile)
	print p.resDict[tvar]
	rtmap = m.getResiTargetMap(p, target)
	if len(rtmap) < 1:
		return
	print 'map %s: %s' % (tvar, repr(rtmap[tvar]))
	return (tvar, rtmap[tvar][0], rtmap[tvar][1])



# function for parsing sdii result
# msai -> seqi -> 'B529(V)'
# pdbseqDict: 132 : 'B529(V)'
# sdiline: [1042-2032-3128 0.006242240179705]
# discard for implant alignment
'''
def sdiiparse(sdiiline, msai2seqi, pdbseqDict):
	split1 = sdiiline.split(' ')
	v_dep = split1[1].strip()

	split2 = split1[0].split('-')
	# some of the indices won't be in the msai2seqi since the column is significant but the position on target pdb msa seq are gaps
	return '%s %s' % ('-'.join([pdbseqDict[msai2seqi[int(msai)] if int(msai) in msai2seqi else -1] for msai in split2]), v_dep)
'''

# convert
# 1042-2032-3128 0.006242240179705
# 1931-2177-3128 0.001309941125401
# 2136-3128-3140 0.003996312858620
# to
#
#
#
# protein.seqDict{} : [132 : 'B529(V)']
# msai -> seqi -> 'B529(V)'
# discard for implant alignment
'''
def sdii2resi():
	if len(sys.argv) < 5:
		print 'sdii2resi: convert msa position to residue number in pdb for a sdii result file' 
		print 'example: python utils_msa.py sdii2resi 1k2p_PF07714_full.fa.3128_3_sdii 1k2p_PF07714_full.fa 1k2p.pdb\n'
		return

	sdiifile = sys.argv[2]
	msafile = sys.argv[3]
	pdbfile = sys.argv[4]

	p = protein(pdbfile)
	m = msa(msafile)
	seqi2msai, msai2seqi = m.getPosMap(p)

	with open(sdiifile) as f:
		sdiilines = f.readlines()

	for line in sdiilines:
		print sdiiparse(line, msai2seqi, p.seqDict)
'''

def writeUniprotSeq():
	if len(sys.argv) < 4:
		print 'writeUniprotSeq: get msa sequence without gaps by searching fasta name'
		print 'example: python utils_msa.py writeuniprotseq PF07714_full.txt BTK_HUMAN\n'
		return

	msafile = sys.argv[2]
	msaheader = sys.argv[3].upper()
	outfile = '%s-%s.seq' % (msafile[0:7], msaheader)


	print 'msa file: %s' % msafile
	print 'target entry: %s' % msaheader

	msaseq = ''
	m = msa(msafile)
	m.setTarget(msaheader)

	for s in m.msaArray:
		if msaheader in s[0]:
			msaheader = s[0]
			msaseq = s[1]

	outputSeq = []
	for a in msaseq:
		if a in ['.', '-', '_']:
			continue
		else:
			outputSeq.append(a)

	print msaheader
	print 'writing: %s' % outfile
	fout = open(outfile, 'w')
	fout.write(''.join(outputSeq))
	fout.close()


def getSeqbyName():
	if len(sys.argv) < 4:
		print 'getSeqbyName: get msa sequence without gaps by searching fasta header (uniprot KB)'
		print 'example: python utils_msa.py getseqbyname PF07714_full.fa BTK_HUMAN\n'
		return

	msafile = sys.argv[2]
	msaheader = sys.argv[3].upper()
	#print 'msa file: %s' % msafile
	#print 'target entry: %s' % msaheader

	msaseq = ''
	m = msa(msafile)
	#m.setTarget(msaheader)

	for s in m.msaArray:
		if msaheader in s[0]:
			msaheader = s[0]
			msaseq = s[1]

	outputSeq = []
	for a in msaseq:
		if a in ['.', '-', '_']:
			continue
		else:
			outputSeq.append(a.upper())

	#print msaheader
	print ''.join(outputSeq)


def getMsabyName():
	if len(sys.argv) < 4:
		print 'getMsabyName: get msa sequence with gaps by searching fasta name'
		print 'example: python utils_msa.py getmsabyname PF07714_full.fa BTK_HUMAN\n'
		return

	msafile = sys.argv[2]
	msaheader = sys.argv[3].upper()
	#print 'msa file: %s' % msafile
	#print 'target entry: %s' % msaheader

	msaseq = ''
	m = msa(msafile)
	#m.setTarget(msaheader)

	for s in m.msaArray:
		if msaheader in s[0]:
			#print s[0]
			print s[1]
			break

def pdistDistribution():
	if len(sys.argv) < 3:
		print 'pdist: write pairwise sequence simiarity value in a file' 
		print 'example: python utils_msa.py pdist PF07714_full.fa BTK_HUMAN\n'
		return

	msafile = sys.argv[2]
	target = sys.argv[3]
	outfile = '%s.pdist' % msafile

	print 'msa file: %s\ntarget header: %s\noutfile: %s\n' % (msafile, target, outfile)

	print 'loading msa ...'
	m = msa(msafile)
	m.setTarget(target)
	print

	print 'saving to [%s] ...' % outfile
	pdist = m.getpdist()
	np.savetxt(outfile, pdist, fmt='%.8', delimiter='\n')
	print '%d pdist saved.' % len(pdist)


def reduceByHamming():
	if len(sys.argv) < 3:
		print 'reduceByHamming: reduce a msa file by selecting sequences that have sequential similarity < 62% (hamming dist > 0.38)'
		print 'example: python utils_msa.py reducebyhamming PF07714_full.fa BTK_HUMAN\n'
		return

	msafile = sys.argv[2]
	target = sys.argv[3]
	outfile = '%s.s62' % msafile

	print 'msa file: %s\ntarget header: %s\noutfile: %s\n' % (msafile, target, outfile)

	print 'loading msa ...'
	m = msa(msafile)
	m.setTarget(target)
	print 

	m.hammingReduction(outfile, 0.62)



def reduceByWeight():
	if len(sys.argv) < 5:
		print 'reduceByWeight: reduce a msa file by weighing and reduce scale (x%)'
		print 'example: python utils_msa.py reducebyweight 1k2p_PF07714_full.fa test.weight BTK_HUMAN 0.5\n'
		return

	msafile = sys.argv[2]
	weightfile = sys.argv[3]
	target = sys.argv[4]
	scale = float(sys.argv[5])
	outfile ='%s.r%d' % (msafile, scale*100)
	print 'msa file: %s' % msafile
	print 'weight file: %s' % weightfile
	print 'target: %s' % target
	print 'reduce scale: %f' % scale
	print 'output file: %s' % outfile

	weight = np.loadtxt(weightfile, delimiter=',')
	print 'weight loaded : %s' % repr(weight.shape)

	print 'loading msa file ...'
	m = msa(msafile)
	m.setTarget(target)

	rlist=[]
	for i in xrange(0, len(weight)):
		rlist.append((i, weight[i]))

	# 0 -> len(weight)
	# small -> large
	sort_rlist = sorted(rlist, key=lambda x: x[1])

	#for k in xrange(0, len(sort_rlist)):
	#	print '[%d]:[%s]' % (k, repr(sort_rlist[k]))

	goal = int(len(weight) * (1-scale))

	target_flag = False
	fout = open(outfile, 'w')
	# save msa sequences with large weights
	print 'Writing output ...'
	for k in xrange(goal, len(weight)):
		(index, w) = sort_rlist[k]
		#print '%d, %f' % (index, w)
		if m.msaArray[index][0] == m.target[0]:
			target_flag = True
		fout.write('>%s\n%s\n' % (m.msaArray[index][0], m.msaArray[index][1]))
	if target_flag == False:
		print 'Inserting target sequence: %s' % m.target[0]
		fout.write('>%s\n%s\n' % (m.target[0], m.target[1]))
	fout.close()
	print 'reduced msa: [%s]\nlen: %d' % (outfile, goal)


def findfamiliar():
	if len(sys.argv)<4:
		print 'findfamiliar: find sequences with mutually hamming distance'
		print 'example: python utils_msa.py findfamiliar PF07714_full.txt BTK_HUMAN 0.8 0.1,0.4'
		return

	msafile = sys.argv[2]
	target = sys.argv[3]
	gap_cutoff = float(sys.argv[4]) # gap cutoff
	range_str = sys.argv[5] # hamming cutoff
	rangeArray = [float(i) for i in range_str.split(',')]
	rangeArray.sort()

	print '\nmsa file: %s' % msafile
	print 'target: %s' % target
	print 'gap cutoff: %f' % gap_cutoff
	print 'hamming cutoff range: %s' % repr(rangeArray)

	print 'loading msa file ...'
	m = msa(msafile)
	m.setTarget(target)

	print 'filtering sequences ...'
	row_set = m.find_familiar(gap_cutoff, rangeArray)

	for i in row_set:
		(head, msaline) = m.msaArray[i]
		seq = msaline.replace('.','').replace('-','').lower()
		headArray = head.split('/')
		outname = '%s-%s-%d.sseq' % (msafile[0:7], headArray[0], i)
		#print '%s %s' % (outname, seq)
		fout = open(outname, 'w')
		fout.write(seq)
		fout.close()

	print 'findsimilar: %d sequences generated.' % len(row_set)


def extractnseq():
	if len(sys.argv)<3:
		print 'findfamiliar: extract n sequences from MSA file' 
		print 'example: python utils_msa.py extractnseq PF07714_full.txt number_of_sequence'
		return

	msafile = sys.argv[2]
	nseq = int(sys.argv[3]) 

	print '\nmsa file: %s' % msafile
	print 'nseq: %d' % nseq

	print 'loading msa file ...'
	m = msa(msafile)
	#m.setTarget(target)

	if len(m.msaArray) < nseq:
		print 'not enough sequence in MSA'
		exit(-1)

	for i in xrange(0,nseq):
		(head, msaline) = m.msaArray[i]
		headArray = head.split('/')
		seq = msaline.replace('.','').replace('-','').lower()
		outname = '%s-%s-%d.sseq' % (msafile[0:7], headArray[0], i)
		#print '%s %s' % (outname, seq)
		fout = open(outname, 'w')
		fout.write(seq)
		fout.close()

	print 'extractnseq: %d sequences generated.' % nseq


def MSAReduction():
	if len(sys.argv) < 4:
		print 'msareduction: reduce columns and rows by cutoffs'
		print 'example: python utils_msa.py msareduction PF07714_full.fa BTK_HUMAN 0.8 0.38'
		print 'python utils_msa.py msareduction PF13855_full.txt RTN4R_HUMAN 0.8 -1'
		return

	msafile = sys.argv[2]
	target = sys.argv[3]
	gap_cutoff = float(sys.argv[4]) # gap cutoff
	hamming_cutoff = float(sys.argv[5]) # hamming cutoff

	print '\nmsa file: %s' % msafile
	print 'target: %s' % target
	print 'gap cutoff: %f' % gap_cutoff
	print 'hamming cutoff: %s' % hamming_cutoff

	print 'loading msa file ...'
	m = msa(msafile)
	m.setTarget(target)

	(seqboard, scoreboard, column_index, row_index) = m.get_msaboard_RC_RR(gap_cutoff, hamming_cutoff)

	'''
	print 'score matrix:'
	for i in xrange(0, len(score)):
		print repr(score[i])
	print 'column index: %s' % repr(column_index)
	print 'row index: %s' % repr(row_index)
	'''
	#seqs = np.array(seqboard)[row_index,:][:,column_index]
	seqs = np.array(seqboard)[row_index,:] # complete column reduced row

	fout = open(msafile+'.rseq', 'w')
	for i in xrange(0, len(row_index)):
		header = m.msaArray[row_index[i]][0]
		fout.write('>'+header+'\n')
		fout.write(''.join(seqs[i,:])+'\n')
	fout.close()
	print 'save reduced sequences to file: [%s]' % (msafile+'.rseq')
	#np.savetxt(msafile+'.rseq', seqs, delimiter='')

	score = np.array(scoreboard)[row_index,:][:,column_index]
	#np.savetxt(msafile+'.score', score, fmt='%.8', delimiter=',')
	print 'save score to file: [%s]' % (msafile+'.score')
	np.savetxt(msafile+'.score', score, fmt='%d', delimiter=',')
	print 'save reduced row indices to file: [%s]' % (msafile+'.row')
	fout = open(msafile+'.row', 'w')
	fout.write(','.join([str(i) for i in row_index])+'\n')
	fout.close()

	print 'save reduced column indices to file: [%s]' % (msafile+'.col')
	fout = open(msafile+'.col', 'w')
	fout.write(','.join([str(i) for i in column_index])+'\n')
	fout.close()



def searchpdbseq():
	if len(sys.argv) < 2:
		print 'searchpdbseq: locate pdb sequence in MSA' 
		print 'example: python utils_msa.py searchpdbseq PF07714_full.fa 1T49_A.pdb\n'
		return	

	msafile = sys.argv[2]
	target = sys.argv[3]

	print 'msa file: %s' % msafile
	print 'pdb target: %s' % target

	m = msa(msafile)
	p = protein(target)

	if m.searchTargetPDB(p)==0:
		print 'cannot locate pdb sequence in MSA'


def resi2msai():
	if len(sys.argv) < 5:
		print 'resi2target: given a residue number output the corresponding position in target msa'
		print 'python utils_msa.py resi2msai PF00014_full.txt BPT1_BOVIN 5pti_pf.pdb A6'
		return

	msafile = sys.argv[2]
	target = sys.argv[3]
	pdbfile = sys.argv[4]
	tvar = sys.argv[5]

	print 'msafile: %s\ntarget header: %s\npdbfile: %s\ntarget variable: %s' % (msafile, target, pdbfile, tvar)
	m = msa(msafile)
	p = protein(pdbfile)
	print p.resDict[tvar]
	rtmap = m.getResiTargetMap(p, target)
	if len(rtmap) < 1:
		return
	print 'map %s: %s' % (tvar, repr(rtmap[tvar]))
	return (tvar, rtmap[tvar][0], rtmap[tvar][1])



# read hcg file

# 5pti_pf.tip,E,A7
# 5pti_pf.tip,P,A8
# 5pti_pf.tip,SD,A47 A50
# 5pti_pf.tip,YA,A21 A48
# 5pti_pf.tip,AIG,A16 A18 A37
# 5pti_pf.tip,QNA,A31 A24 A27

# sequences:
# ..................................FCLE.PP...Y....TG......P.....C.K.......A....................RI........IRYFYN............AKA...G....L.C...QT..........F..V..Y..G.....G.C...R....A...K...R...............NN.F....KSAE..DCMRTC.g.....
# FCLEPPYTGPCKARIIRYFYNAKAGLCQTFVYGGCRAKRNNFKSAEDCMRTCG

# return sorted msa column groups
# [[210, 215], [106, 211], [73, 95, 166], [109, 124, 139]]
def hcg2msa(filename, rtmap):
	fin = open(filename, 'r')
	lines = fin.readlines()
	fin.close()

	msaGroupArray = []
	for line in lines:
		line = line.strip()
		if len(line)<1:
			continue
		print line
		lineArray = line.split(',')
		if len(lineArray[1])<2:
			continue
		# remove chain ID
		#resiStr = filter(lambda x: x.isdigit() or x==' ', lineArray[2])
		resiStr = lineArray[2] # A31 A24 A27
		#print resiStr.split(' ')
		msaColArray = [rtmap[resi][0] for resi in resiStr.split(' ')] # [139, 109, 124]
		#print msaColArray
		msaColArray.sort()
		msaGroupArray.append(msaColArray)
		#msaGroupArray.append('-'.join([str(i) for i in msaColArray])) # ['210-215', '106-211', '73-95-166', '109-124-139']

		#print msaColArray
		#print ''

	print repr(msaGroupArray)
	return msaGroupArray


# load sdii into a dictionary
'''
210-215 0.318813279927245
106-211 0.151691013032345
73-95-166 0.056163522315670
109-124-139 0.175665187777980
65-135-139 0.000047065100507
135-151 0.009274777159306
'''
# return dictionary with (key: 39-140-210, value = 0.0788593466276019)
def loadsdii(sdiifile):
	fin = open(sdiifile, 'r')
	lines = fin.readlines()
	fin.close()

	sdict = {}
	for line in lines:
		line = line.strip()
		if len(line)<1:
			continue
		strArray = line.split(' ')
		sdict[strArray[0]] = float(strArray[1])
	
	#print repr(sdict)
	return sdict	

# calculate columwised substritution marix
# input a column of MSA in list type
def calcColSM(sm, msaMatrix, column):
	#rc = msaMatrix[:, column].replace('-','.').replace('O','.').replace('U','.').replace('X','.').replace('B','.').replace('Z','.')
	colset = set(msaMatrix[:, column])
	if '.' in colset:
		colset.remove('.')
	if '-' in colset:
		colset.remove('-')
	if 'O' in colset:
		colset.remove('O')
	if 'U' in colset:
		colset.remove('U')
	if 'X' in colset:
		colset.remove('X')
	if 'B' in colset:
		colset.remove('B')
	if 'Z' in colset:
		colset.remove('Z')
	#alphabet = sorted(set(msaMatrix[:, column])) # get sorted alphabet list
	alphabet = sorted(colset) # get sorted alphabet list
	#print repr(alphabet)
	# get A frequency
	freqDict = collections.Counter(msaMatrix[:, column]) # get frequency
	#print repr(freqDict)

	# calculate AA frequency cij
	for i in xrange(0, len(alphabet)):
		A = alphabet[i].upper()
		if A not in alphabet:
			continue
		for j in xrange(i+1, len(alphabet)):
			B = alphabet[j].upper()
			if B not in alphabet:
				continue
			sm[A+B] += freqDict[A]*freqDict[B]
		sm[A+A] += freqDict[A]*(freqDict[A]-1)/2

	return 


# read hcg pdb msa(raw) sdii pdbtitle
# output new substitute matrix
def hcg2blossum():
	if len(sys.argv) < 4:
		print 'hcg2blossum: construct new substitution matrix from contact group'
		print 'example:python utils_msa.py hcg2blossum 5pti_pf.pdb 5pti_pf.tip.hcg PF00014_full.txt.rseq PF00014_full.txt.sdii PF00014_full.txt.sm BPT1_BOVIN'
		print 'output: a substitution matrix file (same format as BLOSSUM62)'
		return
	#print sys.argv[0] # utils_msa.py
	#print sys.argv[1] # hcg2blossum
	pdbfile = sys.argv[2] # pdb name
	hcgfile = sys.argv[3] # hcg
	msafile = sys.argv[4] # msa (full or reduced)
	sdiifile = sys.argv[5] # sdii
	outfile = sys.argv[6] # new substitution matrix
	target = sys.argv[7] # target name

	# get msa in matrix format
	m = msa(msafile)
	msaMatrix = np.array([list(s[1]) for s in m.msaArray]) # matrix format of msa

	#for i in xrange(0, len(seqs)):
	#	print seqs[i]
	print 'msa matrix: ' + repr(msaMatrix.shape)

	# get resi -> msai map	
	p = protein(pdbfile)
	rtmap = m.getResiTargetMap(p, target)

	sdiidict = loadsdii(sdiifile) # key: 39-140-210, value = 0.0788593466276019
	msaGroupArray = hcg2msa(hcgfile, rtmap) # [[210, 215], [106, 211], [73, 95, 166], [109, 124, 139]]

	# init substitution matrix
	EBlist = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X', '*']
	#AAlist = sorted(EBlist)
	#AAlist = sorted(['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X', '*'])
	AAlist = sorted(['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'])
	sm = {}
	for i in xrange(0, len(AAlist)):
		for j in xrange(i, len(AAlist)):
			key = '%s%s' % (AAlist[i], AAlist[j])
			sm[key] = 0
	print AAlist
	print 'Alphabet: %d' % len(AAlist) 
	print 'AA: %d' % len(sm)

	# accumulate substitution matrix AA frequency for all the contact group columns
	# Sum the scores for each columns across column
	print ''
	w = 0 # count column number
	for mg in msaGroupArray:
		# form key for co-evolve value 
		sdiikey = '-'.join([str(i) for i in mg])
		if sdiikey not in sdiidict:
			print 'hcg2blossum():discard group: %s' % sdiikey
			continue
		sdiiweight = sdiidict[sdiikey]
		print (sdiikey, sdiiweight)

		# accumulate SM for each contact group / column group
		for col in mg:
			w +=1
			calcColSM(sm, msaMatrix, col)
			# debug
			'''
			print 'col: %d' % col
			if count==2:
				print repr(sm)
				return
			'''
		print ''

	#print repr(sm)
	#print ''

	n = msaMatrix.shape[0]
	T = w*n*(n-1)/2 # normalization term
	print 'w: %d' % w # number of columns (contact group)
	print 'n: %d' % n # number of sequence
	print 'T: %d' % T


	# convert cij to qij
	# Normalize the pair frequencies so they will sum to 1
	for c in sm:
		sm[c] = 1.0*sm[c]/T

	#print repr(sm)
	#print ''

	# Calculate the expected probability of occurrence of the ith residue in an (i,j) pair
	# pi = qii + sum( qij/2 )_{i!=j}
	pi = {}
	for i in xrange(0, len(AAlist)):
		A = AAlist[i]
		sum_qij = 0
		for j in xrange(i+1, len(AAlist)): # i should not = j
			B = AAlist[j]
			sum_qij += sm[A+B]/2
		pi[A] = sm[A+A] + sum_qij

	print repr(pi)	
	print ''

	# The desired denominator is the expected frequency for each pair 
	eij = {}
	for i in xrange(0, len(AAlist)):
		A = AAlist[i]
		for j in xrange(i+1, len(AAlist)):
			B = AAlist[j]
			eij[A+B] = 2 * pi[A] * pi[B]
		eij[A+A] = pi[A] * pi[A]

	print len(eij)
	print repr(eij)	
	print ''

	#  Log odds ratio sij = round(2*log2(qij/eij))
	sij = {}
	for i in xrange(0, len(AAlist)):
		A = AAlist[i]
		for j in xrange(i, len(AAlist)):
			B = AAlist[j]
			if eij[A+B] == 0.0 or sm[A+B]==0.0:
				sij[A+B] = 0
			else:
				sij[A+B] = int(round(2*math.log((sm[A+B]/eij[A+B]),2)))
			#	sij[A+B] = sm[A+B]/eij[A+B]

	print repr(sij)	
	print len(sij)
	print ''

	saveBlosum(EBlist, sij, outfile)


# write new substitution matrix to a file
# alphabet = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X', '*']
# s: dictionary for a corresponding alphabet
# filename: output filename
def saveBlosum(alphabet, s, filename):
	print 'writing: %s' % filename
	fout = open(filename, 'w')
	for i in xrange(0, len(alphabet)):
		A = alphabet[i]
		lineArray = []
		for j in xrange(0, len(alphabet)):
			B = alphabet[j]
			if A+B in s:
				lineArray.append(str(s[A+B]).rjust(2))
			elif B+A in s:
				lineArray.append(str(s[B+A]).rjust(2))
			else:
				lineArray.append(' 0')
		fout.write(' '.join(lineArray)+'\n')
	fout.close()


# called by ncg2blossum
# convert ncg to sorted msa group array
# example of ncg file:
'''
A93 A36 A89 A37 A35
A97 A23 A96 A24 A27
'''
# rtmap['B641'] = (seqpos, 'R')
def ncg2msa(filename, rtmap):
	fin = open(filename, 'r')
	lines = fin.readlines()
	fin.close()

	#for k in rtmap:
	#	print (k, rtmap[k])
	msaGroupArray = []
	for line in lines:
		line = line.strip()
		if len(line)<1:
			continue
		#print line
		resiStr = line

		msaColArray = []
		discard_flag=False # for inconsistent pfam index. v29.0 1APS 1-97 but in v30.0 1APS 10-97
		for resi in resiStr.split(' '):
			if resi not in rtmap:
				discard_flag = True # if pdb residue not exist in MSA then discard the current contact group mapping
				break
			else:
				msaColArray.append(rtmap[resi][0])

		if discard_flag == True:
			continue
		else:
			msaGroupArray.append(msaColArray)

		#msaColArray = [rtmap[resi][0] for resi in resiStr.split(' ')] # [139, 109, 124]
		#msaColArray.sort() #Do not sort!

	#print repr(msaGroupArray)
	return msaGroupArray


# read hcg pdb msa(raw) sdii pdbtitle
# output new substitute matrix
def ncg2blossum():
	if len(sys.argv) < 7:
		print 'ncg2blossum: construct new substitution matrix from contact group'
		print 'example:python utils_msa.py ncg2blossum 5pti_pf.pdb 5pti_pf.tip.ncg PF00014_full.txt.rseq PF00014_full.txt.sdii BPT1_BOVIN order'
		print 'output: a substitution matrix file (same format as BLOSSUM62)'
		return
	#print sys.argv[0] # utils_msa.py
	#print sys.argv[1] # hcg2blossum
	pdbfile = sys.argv[2] # pdb name
	ncgfile = sys.argv[3] # hcg
	msafile = sys.argv[4] # msa (full or reduced)
	sdiifile = sys.argv[5] # sdii
	target = sys.argv[6] # target name
	order = int(sys.argv[7])
	outfile = msafile[0:7]+".sm" # new substitution matrix

	# get msa in matrix format
	m = msa(msafile)
	msaMatrix = np.array([list(s[1]) for s in m.msaArray]) # matrix format of msa

	#for i in xrange(0, len(seqs)):
	#	print seqs[i]
	print 'msa matrix: ' + repr(msaMatrix.shape)

	# get resi -> msai map	
	p = protein(pdbfile)
	rtmap = m.getResiTargetMap(p, target)

	sdiidict = loadsdii(sdiifile) # key: 39-140-210, value = 0.0788593466276019
	msaGroupArray = ncg2msa(ncgfile, rtmap) # [[210, 215], [106, 211], [73, 95, 166], [109, 124, 139]]

	# get non overlapped column indices
	colset = set()
	for g in msaGroupArray:
		rg = g[0:order] # get ith order contact group
		rg.sort() # for generating key
		sdiikey = '-'.join([str(r) for r in rg])
		if sdiikey not in sdiidict:
			#print 'ncg2sdiicol(): discard group: %s for low sdii' % sdiikey
			continue
		print (sdiikey, sdiidict[sdiikey])			
		for resi in rg: # for significant ncg, add corresponding MSA column index
			colset.add(resi)

	# init substitution matrix
	EBlist = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X', '*']
	#AAlist = sorted(EBlist)
	#AAlist = sorted(['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X', '*'])
	AAlist = sorted(['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'])
	sm = {}
	for i in xrange(0, len(AAlist)):
		for j in xrange(i, len(AAlist)):
			key = '%s%s' % (AAlist[i], AAlist[j])
			sm[key] = 0
	print AAlist
	print 'Alphabet: %d' % len(AAlist) 
	print 'AA: %d' % len(sm)

	# accumulate substitution matrix AA frequency for all the contact group columns
	# Sum the scores for each columns across column
	print ''
	w = 0 # count column number
	for col in colset:
		w+=1
		calcColSM(sm, msaMatrix, col)
	'''
	for mg in msaGroupArray:
		# form key for co-evolve value 
		sdiikey = '-'.join([str(i) for i in mg])
		if sdiikey not in sdiidict:
			print 'hcg2blossum():discard group: %s' % sdiikey
			continue
		sdiiweight = sdiidict[sdiikey]
		print (sdiikey, sdiiweight)

		# accumulate SM for each contact group / column group
		for col in mg:
			w +=1
			calcColSM(sm, msaMatrix, col)
		print ''
	'''
	#print repr(sm)
	#print ''

	n = msaMatrix.shape[0]
	T = w*n*(n-1)/2 # normalization term
	print 'w: %d' % w # number of columns (contact group)
	print 'n: %d' % n # number of sequence
	print 'T: %d' % T


	# convert cij to qij
	# Normalize the pair frequencies so they will sum to 1
	for c in sm:
		sm[c] = 1.0*sm[c]/T

	#print repr(sm)
	#print ''

	# Calculate the expected probability of occurrence of the ith residue in an (i,j) pair
	# pi = qii + sum( qij/2 )_{i!=j}
	pi = {}
	for i in xrange(0, len(AAlist)):
		A = AAlist[i]
		sum_qij = 0
		for j in xrange(i+1, len(AAlist)): # i should not = j
			B = AAlist[j]
			sum_qij += sm[A+B]/2
		pi[A] = sm[A+A] + sum_qij

	print repr(pi)	
	print ''

	# The desired denominator is the expected frequency for each pair 
	eij = {}
	for i in xrange(0, len(AAlist)):
		A = AAlist[i]
		for j in xrange(i+1, len(AAlist)):
			B = AAlist[j]
			eij[A+B] = 2 * pi[A] * pi[B]
		eij[A+A] = pi[A] * pi[A]

	print len(eij)
	print repr(eij)	
	print ''

	#  Log odds ratio sij = round(2*log2(qij/eij))
	sij = {}
	for i in xrange(0, len(AAlist)):
		A = AAlist[i]
		for j in xrange(i, len(AAlist)):
			B = AAlist[j]
			if eij[A+B] == 0.0 or sm[A+B]==0.0:
				sij[A+B] = 0
			else:
				sij[A+B] = int(round(2*math.log((sm[A+B]/eij[A+B]),2)))
			#	sij[A+B] = sm[A+B]/eij[A+B]

	print repr(sij)	
	print len(sij)
	print ''

	saveBlosum(EBlist, sij, outfile)


# called in ncg2sdiicol()
# load map between pdb residue ID and MSA uniprot position ID 
# dictionary element: ('A9', (14, 'V')) : (chain+resi, (MSA position index, resn))
# mapfile:
# AT284 1218 T  : chain A residue T resn 284 => position 1218 (start from 0) resn T
# AE285 1220 e  : lowercase exists!
# AR286 -1 -	: residue number that cannot map to MSA position (does not exist)
def getPDBUniprotMap(mapfile):
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


# by just using top sdii columns alone to save sdiicol file
# no structural info considered
def topsdii2sdiicol():
	if len(sys.argv)<3:
		print 'topsdii2sdiicol: selecting MSA column with top sdii value only into .sdiicol file'
		print 'python utils_msa.py topsdii2sdiicol PF00589_full.txt.rseq PF00589_full.txt.all_2_sdii.4.top'
		return	

	msafile = sys.argv[2]
	sdiifile = sys.argv[3]
	outfile = msafile[0:7]+'.tsdiicol'

	print 'msafile :%s' % msafile
	print 'sdiifile :%s' % sdiifile

	colset = set()
	with open(sdiifile) as f:
		for line in f:
			if len(line)<2:
				continue
			sdiiArray = line.strip().split(' ')
			for i in sdiiArray[0].split('-'):
				colset.add(int(i))
	print 'topsdii2sdiicol():writing %s, sdiicol %d' % (outfile, len(colset))
	print '%s' % (repr(sorted(colset)))

	fout = open(outfile, 'w')
	fout.write('%s %s\n' % (msafile, ' '.join([str(c) for c in sorted(colset)])))
	fout.close()

#by just using residue contact to save msa columns
def cg2col():
	if len(sys.argv)<5:
		print 'cg2col: write selected MSA column into .csdiicol file'
		print 'python utils_msa.py cg2col 1a0p-A.rpdb.tip.ncg PF00589_full.txt.rseq 1a0p-A-PF00589-XERD_ECOLI.map 2'
		return

	ncgfile = sys.argv[2] # hcg
	msafile = sys.argv[3] # msa (full or reduced)
	resimapfile = sys.argv[4]
	orderlist = [int(i) for i in sys.argv[6].split(',')]
	outfile =  msafile[0:7]+'.csdiicol'


	print 'ncgfile :%s' % ncgfile
	print 'msafile :%s' % msafile
	print 'resimapfile :%s' % resimapfile
	print 'ncg order list : [%s]' % repr(orderlist)

	# get msa in matrix format
	m = msa(msafile)
	msaMatrix = np.array([list(s[1]) for s in m.msaArray]) # matrix format of msa

	#for i in xrange(0, len(seqs)):
	#	print seqs[i]
	print 'msa matrix: ' + repr(msaMatrix.shape)

	#rtmap = m.getResiTargetMap(p, target) # ('A9', (14, 'V')) : (resi+chain, (MSA index, resn))
	rtmap = getPDBUniprotMap(resimapfile) # ('A9', (14, 'V')) : (resi+chain, (MSA index, resn))

	sdiidict = loadsdii(sdiifile) # key: 39-140-210, value = 0.0788593466276019
	msaGroupArray = ncg2msa(ncgfile, rtmap) # unsorted [[86, 83, 198, 127, 120], [138, 76, 82, 127, 132]]

	# output msa column set
	colset = set()
	for i in orderlist:
		for g in msaGroupArray:
			rg = g[0:i] # get ith order contact group
			rg.sort() # for generating key
			#sdiikey = '-'.join([str(r) for r in rg])
			#if sdiikey not in sdiidict:
				#print 'ncg2sdiicol(): discard group: %s for low sdii' % sdiikey
			#	continue
			#print (sdiikey, sdiidict[sdiikey])			
			for resi in rg: # for significant ncg, add corresponding MSA column index
				colset.add(resi)

	print 'cg2col():writing %s, col %d' % (outfile, len(colset))
	print '%s' % (repr(sorted(colset)))
	fout = open(outfile, 'w')
	#fout.write(' '.join([str(c) for c in colset]))
	fout.write('%s %s\n' % (msafile, ' '.join([str(c) for c in colset])))
	fout.close()



# from possible residue contact
# select informative sdii columns from MSA and save the column indices into file
# output: .sdiicol file
def ncg2sdiicol():
	if len(sys.argv)<6:
		print 'ncg2sdiicol: write selected MSA column into .sdiicol file'
		print '$ python utils_msa.py ncg2sdiicol 1a0p-A.rpdb.tip.ncg PF00589_full.txt.rseq PF00589_full.txt.all_2_sdii.4.top 1a0p-A-PF00589-XERD_ECOLI.map 2'
		return

	ncgfile = sys.argv[2] # hcg
	msafile = sys.argv[3] # msa (full or reduced)
	sdiifile = sys.argv[4] # sdii
	resimapfile = sys.argv[5]
	orderlist = [int(i) for i in sys.argv[6].split(',')]
	outStringArray = sdiifile.split('.')
	outfile =  msafile[0:7]+'.'+outStringArray[-1]+'.msacol' # new substitution matrix


	print 'ncgfile :%s' % ncgfile
	print 'msafile :%s' % msafile
	print 'sdiifile :%s' % sdiifile
	print 'resimapfile :%s' % resimapfile
	print 'ncg order list : [%s]' % repr(orderlist)

	# get msa in matrix format
	m = msa(msafile)
	msaMatrix = np.array([list(s[1]) for s in m.msaArray]) # matrix format of msa

	#for i in xrange(0, len(seqs)):
	#	print seqs[i]
	print 'msa matrix: ' + repr(msaMatrix.shape)

	#rtmap = m.getResiTargetMap(p, target) # ('A9', (14, 'V')) : (resi+chain, (MSA index, resn))
	rtmap = getPDBUniprotMap(resimapfile) # ('A9', (14, 'V')) : (resi+chain, (MSA index, resn))

	sdiidict = loadsdii(sdiifile) # key: 39-140-210, value = 0.0788593466276019
	msaGroupArray = ncg2msa(ncgfile, rtmap) # unsorted [[86, 83, 198, 127, 120], [138, 76, 82, 127, 132]]

	# output msa column set
	colset = set()
	for i in orderlist:
		for g in msaGroupArray:
			rg = g[0:i] # get ith order contact group
			rg.sort() # for generating key
			sdiikey = '-'.join([str(r) for r in rg])
			if sdiikey not in sdiidict:
				#print 'ncg2sdiicol(): discard group: %s for low sdii' % sdiikey
				continue
			print (sdiikey, sdiidict[sdiikey])			
			for resi in rg: # for significant ncg, add corresponding MSA column index
				colset.add(resi)

	print 'ncg2sdiicol():writing %s, col %d' % (outfile, len(colset))
	print '%s' % (repr(sorted(colset)))
	fout = open(outfile, 'w')
	#fout.write(' '.join([str(c) for c in colset]))
	fout.write('%s %s\n' % (msafile, ' '.join([str(c) for c in colset])))
	fout.close()


# read combined .sdiicol file
# for each one calculate the substitution rate cij
def sdii2blosum():
	if len(sys.argv) < 2:
		print 'sdii2blosum: construct new substitution matrix from sdii filtered msa columns'
		print 'python utils_msa.py sdii2blosum PF00014.sdiicol'
		print 'output: a substitution matrix file (same format as BLOSSUM62)'
		return

	sdiicolfile = sys.argv[2]
	outfile = sdiicolfile+".sm" # new substitution matrix

	sdiicolArray = []
	with open(sdiicolfile) as fp:
		for line in fp:
			line = line.strip()
			if len(line) < 1:
				print 'error sdiicol line: %s' % line

			valueArray = line.split(' ')
			msafile = valueArray[0]
			if len(valueArray)==1:
				print '%s: no column' % msafile
				continue
			colset = [int(i) for i in valueArray[1:]]
			sdiicolArray.append((msafile, colset))
	print 'sdii2bolsum(): %d sdiicol lines loaded ..' % len(sdiicolArray)

	# init substitution matrix
	AAlist = sorted(['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'])
	sm = {}
	for i in xrange(0, len(AAlist)):
		for j in xrange(i, len(AAlist)):
			key = '%s%s' % (AAlist[i], AAlist[j])
			sm[key] = 0
	#print AAlist
	#print 'Alphabet: %d' % len(AAlist) 
	#print 'AA: %d' % len(sm)

	T=0 # total number of possible substitutions
	for msafile, colset in sdiicolArray:
		#print len(colset), repr(colset)

		# get msa in matrix format
		m = msa(msafile)
		msaMatrix = np.array([list(s[1]) for s in m.msaArray]) # matrix format of msa
		#print 'msa matrix: ' + repr(msaMatrix.shape)

		# accumulate substitution matrix AA frequency for all the contact group columns
		# Sum the scores for each columns across column
		for col in colset:
			calcColSM(sm, msaMatrix, col)

		w = len(colset) # number of columns for current MSA
		n = msaMatrix.shape[0]
		T= sum(sm.itervalues())
		#T = w*n*(n-1)/2 # normalization term, removing gaps will cause T in action is smaller than w*n*(n-1)/2
		print '%s: w: %d/%d, n: %d, T: %d' % (msafile, w, msaMatrix.shape[1], n, T)
		#T+=T1

	print 'total T:  %d' % T
	print 'total sm: %d' % sum(sm.itervalues())
	print ''
	# convert cij to qij
	# Normalize the pair frequencies so they will sum to 1
	for c in sm:
		sm[c] = 1.0*sm[c]/T
	print 'sum(sm): %f' % sum(sm.itervalues())
	print repr(sm)[0:10]
	print ''

	# Calculate the expected probability of occurrence of the ith residue in an (i,j) pair
	# pi = qii + sum( qij/2 )_{i!=j}
	pi = {}
	for i in xrange(0, len(AAlist)):
		A = AAlist[i]
		sum_qij = 0
		for j in xrange(i+1, len(AAlist)): # i should not = j
			B = AAlist[j]
			sum_qij += sm[A+B]/2
		pi[A] = sm[A+A] + sum_qij
	print 'sum(pi): %f' % sum(pi.itervalues())
	#print repr(pi)[0:10]	
	print ''

	# The desired denominator is the expected frequency for each pair 
	eij = {}
	for i in xrange(0, len(AAlist)):
		A = AAlist[i]
		for j in xrange(i+1, len(AAlist)):
			B = AAlist[j]
			eij[A+B] = 2 * pi[A] * pi[B]
		eij[A+A] = pi[A] * pi[A]

	print 'sum(eij): %f' % sum(eij.itervalues())
	#print repr(eij)[0:10]
	print ''

	#  Log odds ratio sij = round(2*log2(qij/eij))
	sij = {}
	for i in xrange(0, len(AAlist)):
		A = AAlist[i]
		for j in xrange(i, len(AAlist)):
			B = AAlist[j]
			if eij[A+B] == 0.0 or sm[A+B]==0.0:
				sij[A+B] = 0
			else:
				sij[A+B] = int(round(2*math.log((sm[A+B]/eij[A+B]),2)))
			#	sij[A+B] = sm[A+B]/eij[A+B]

	print repr(sij)	
	print len(sij)
	print ''

	EBlist = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X', '*']
	saveBlosum(EBlist, sij, outfile)



#
def applysm():
	if len(sys.argv) < 6:
		print 'applysm: write new sm file'
		print 'python utils_msa.py applysm untitled_blosum62.txt smlist.txt 0.1 1 outfile'
		return
	alphabet = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X', '*']

	blosumfile = sys.argv[2]
	newsmfilelist = sys.argv[3]
	weight = float(sys.argv[4])
	title = int(sys.argv[5])


	fin = open(newsmfilelist, 'r')
	lines = fin.readlines()
	fin.close()

	newsm = np.loadtxt(lines[0].strip())
	for i in xrange(1, len(lines)):
		line = lines[i].strip()
		if len(line) < 1:
			continue
		print 'reading %s' % line
		newsm += np.loadtxt(line)
	origblosum = np.loadtxt(blosumfile)

	newblosum = np.rint(origblosum + weight * newsm)

	#outfilename = 'combined-'+blosumfile+'.'+sys.argv[4]+'.sm'
	outfilename = sys.argv[6]
	if title==1:
		fout = open(outfilename, 'w')
		fout.write('   '+'  '.join(alphabet)+'\n')
		for i in xrange(0, newblosum.shape[0]):
			fout.write(alphabet[i]+' '+' '.join(['%2i'%n for n in newblosum[i,:]])+'\n')
		fout.close()
		print 'write: %s with title' % (outfilename)
	else:
		np.savetxt(outfilename, newblosum, fmt='%2i', delimiter=' ')
		print 'write: %s without title' % (outfilename)


#######################################################################################################
def main():

	dispatch = {
		'resi2msai': resi2msai, 'msai2resi':msai2resi, 'sdii2resi': sdii2resi, 'getseqbyname': getSeqbyName, 'getmsabyname': getMsabyName,
		'reducebyweight': reduceByWeight, 'reducebyhamming': reduceByHamming, 'resi2target': resi2target, 'pdist': pdistDistribution, 'msareduction':MSAReduction,
		'searchpdbseq': searchpdbseq, 'hcg2blossum': hcg2blossum, 'applysm': applysm, 'ncg2sdiicol':ncg2sdiicol, 'ncg2blossum':ncg2blossum,
		'writeuniprotseq':writeUniprotSeq, 'printuniprot':printUniprot, 'sdii2blosum':sdii2blosum, 'findfamiliar':findfamiliar,'extractnseq':extractnseq,
		'topsdii2sdiicol':topsdii2sdiicol
	}

	if len(sys.argv)<2:
		for k in dispatch:
			dispatch[k]()
		return

	cmd = sys.argv[1]

	flag = False
	for key in dispatch:
		if key == cmd:
			dispatch[key]()
			flag = True
	if flag == False:
		print 'Wrong cmd string'



if __name__ == '__main__':
	main()
