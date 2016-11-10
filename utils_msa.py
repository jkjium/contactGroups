'''
get msa position id 
'''
import sys
import numpy as np
import collections
import math
from msa import msa
from protein import protein

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


def getSeqbyName():
	if len(sys.argv) < 4:
		print 'getSeqbyName: get msa sequence without gaps by searching fasta name'
		print 'example: python utils_msa.py getseqbyname PF07714_full.fa BTK_HUMAN\n'
		return

	msafile = sys.argv[2]
	msaheader = sys.argv[3].upper()
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
	print ''.join(outputSeq)


def getMsabyName():
	if len(sys.argv) < 4:
		print 'getMsabyName: get msa sequence with gaps by searching fasta name'
		print 'example: python utils_msa.py getmsabyname PF07714_full.fa BTK_HUMAN\n'
		return

	msafile = sys.argv[2]
	msaheader = sys.argv[3].upper()
	print 'msa file: %s' % msafile
	print 'target entry: %s' % msaheader

	msaseq = ''
	m = msa(msafile)
	m.setTarget(msaheader)

	for s in m.msaArray:
		if msaheader in s[0]:
			print s[0]
			print s[1]

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


def MSAReduction():
	if len(sys.argv) < 4:
		print 'msareduction: reduce columns and rows by cutoffs'
		print 'example: python utils_msa.py msareduction PF07714_full.fa BTK_HUMAN 0.2 0.62\n'
		return

	msafile = sys.argv[2]
	target = sys.argv[3]
	gap_cutoff = float(sys.argv[4]) # gap cutoff
	hamming_cutoff = float(sys.argv[5]) # hamming cutoff

	print 'msa file: %s' % msafile
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
	colset = set(msaMatrix[:, column])
	if '.' in colset:
		colset.remove('.')
	if 'X' in colset:
		colset.remove('X')
	if 'B' in colset:
		colset.remove('B')
	if 'Z' in colset:
		colset.remove('Z')
	#alphabet = sorted(set(msaMatrix[:, column])) # get sorted alphabet list
	alphabet = sorted(colset) # get sorted alphabet list
	print repr(alphabet)
	# get A frequency
	freqDict = collections.Counter(msaMatrix[:, column]) # get frequency
	print repr(freqDict)

	# calculate AA frequency cij
	for i in xrange(0, len(alphabet)):
		A = alphabet[i]
		for j in xrange(i+1, len(alphabet)):
			B = alphabet[j]
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



#
def applysm():
	if len(sys.argv) < 5:
		print 'applysm: write new sm file'
		print 'python utils_msa.py applysm blosum62.txt smlist.txt 1 1'
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

	outfilename = blosumfile+'.'+sys.argv[4]+'.new'
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
		'searchpdbseq': searchpdbseq, 'hcg2blossum': hcg2blossum, 'applysm': applysm
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
