'''
should stop using this script.
combined with utils_testcase.py
'''

import sys
import subprocess as sp 
import numpy as np
import math
import time
import sys
import commp as cp
from operator import itemgetter

'''
$1:  file name
$2: proram
$3:	 matrix
$4:  gap open penalty
$5:  gap extend penalty
$6:  aligned sequence length
$7:  identity number
$8:  identity percentile 
$9:  similarity number
$10:  similarity percentile
$11:  gaps number
$12: gaps percentile
$13: align score
$14: seq A pure length
$15: aligned seq A
$16: seq B pure length
$17: aligned seq B
'''
class palign(object):
	def __init__(self, flatStr):
		flatArray = flatStr.split(' ')
		if len(flatArray)!= 17:
			print 'Invalid flat string len: %d\n' % len(flatArray)
			print flatStr
			exit(1)
		self.flatstr = flatStr
		self.name = flatArray[0]
		self.program = flatArray[1]
		self.matrix = flatArray[2]
		self.gapopen = float(flatArray[3])
		self.gapextend = float(flatArray[4])
		self.alignlen = int(flatArray[5])
		self.nid = float(flatArray[6])
		self.pid = float(flatArray[7])
		self.nsm = float(flatArray[8])
		self.psm = float(flatArray[9])
		self.ngp = float(flatArray[10])
		self.pgp = float(flatArray[11])
		self.score = float(flatArray[12])
		self.seqAlen = float(flatArray[13])
		self.seqA = flatArray[14]
		self.seqBlen = float(flatArray[15])
		self.seqB = flatArray[16]

	# return the names of the paired sequences
	# p.1aoe.1kmv.seq
	def pairnames(self):
		strArr = self.name.split('.')
		return (strArr[1], strArr[2])

	# return index list for aligned positions
	# non-gap positions
	def alnpos(self):
		gap = ['.', '-']
		poslist = []
		for i in xrange(0, self.alignlen):
			if (self.seqA[i] not in gap) and (self.seqB[i] not in gap):
				poslist.append(i)
		return poslist

	def alnposlist(self):
		gap = ['.', '-']
		setlist = []
		start = 0
		forward = 0
		for i in xrange(0, self.alignlen):
			if (self.seqA[i] not in gap) and (self.seqB[i] not in gap):
				forward = i
			else:
				if forward > start:
					setlist.append([k for k in xrange(start+1, forward+1)])
				start = i 
		#print repr(setlist)
		return setlist


	# dump object to stdout
	def dump(self):
		print '\n-----------------------------'
		print 'name: %s' % self.name
		print 'program: %s' % self.program
		print 'matrix: %s' % self.matrix
		print 'gap open penalty: %f' % self.gapopen
		print 'gap extend penalty: %f' % self.gapextend
		print 'alignment length: %d' % self.alignlen
		print 'number of identity: %d' % self.nid
		print 'percent of identity: %f' % self.pid
		print 'number of similarity: %d' % self.nsm
		print 'percent of similarity: %f' % self.psm
		print 'number of gaps: %d' % self.ngp
		print 'percent of gaps: %f' % self.pgp
		print 'alignment score: %f' % self.score
		print 'seq A length: %d' % self.seqAlen
		print 'aligned seq A: %s' % self.seqA
		print 'seq B length: %d' % self.seqBlen
		print 'aligned seq B: %s' % self.seqB
		print '-----------------------------\n'


class alignflat(object):

	def __init__(self, flatfile):
		self.name = flatfile
		self.flatArray = []
		self.totalnid = 0
		self.totalnsm = 0
		self.totalngp = 0
		self.totalres = 0
		with open(flatfile) as f:
			for line in f:
				if len(line)<2:
					continue
				p = palign(line.strip())
				self.totalnid+=p.nid
				self.totalnsm+=p.nsm
				self.totalngp+=p.ngp
				
				self.totalres+=p.seqAlen
				self.totalres+=p.seqBlen
				#p.dump()
				self.flatArray.append(p)
		#print '%s: %d alignments loaded' % (self.name, len(self.flatArray))

	def dump(self):
		for f in self.flatArray:
			f.dump()

# return number of tp and fp for one blast cath output
def cathtpfp(blastoutfile):
	# cath.3-90-25-10.fa.b62.out
	cathid = blastoutfile.split('.')[1]
	tp = fp = 0
	with open(blastoutfile) as fd:
		for line in fd:
			# 1aq0A00 3-20-20-80,0.0
			line = line.strip()
			if len(line) == 0:
				continue
			s1 = line.split(',')[0]
			outcathid = s1.split(' ')[1]
			if outcathid == cathid:
				tp+=1
			else:
				fp+=1
	if (tp + fp) == 0: # for empty file
		ret = (cathid, (-1,-1))
	else:
		ret = (cathid, (tp, fp))
	return ret


# output a column of tp,fp with stub order
def blastcathbystub(arglist):
	if len(arglist) < 3:
		cp._err('Usage: python utils_testcase.py blastcathbystub stubfile blastoutfilelist outfile')

	stubfile = arglist[0]
	ofilelistfile = arglist[1]
	outfile = arglist[2]

	# load stub into a list
	stub = []
	with open(stubfile) as fp:
		# cath.1-10-10-10.fa
		for line in fp:
			line = line.strip()
			if len(line)==0:
				continue
			sarr = line.split('.')
			stub.append(sarr[1])
	stub.sort() # sort cath family id
	print 

	# load blast output into a dictionary
	ofilelist = []
	with open(ofilelistfile) as fp:
		for line in fp:
			line = line.strip()
			if len(line) == 0:
				continue
			ofilelist.append(line) # PF07582.fa.B62.out

	# for each blast out file calculate tp, fp then store in a dictionary
	# key: pfamid, value: (tp, fp)
	blastout = dict(cathtpfp(file) for file in ofilelist)

	# write outfile with stub order
	with open(outfile, 'w') as fp:
		fp.write('\n'.join(['%s %d %d' % (cathid, blastout[cathid][0], blastout[cathid][1]) for cathid in stub]))
	cp._info('save to %s' % outfile)



# parse blast output ($ blastp -query t.fa -db $PFAM -outfmt "10 stitle evalue" -evalue 0.0001 -matrix BLOSUM80 -o t.fa.out)
# A0A176VM62_MARPO/110-314 A0A176VM62.1 PF02263.18;GBP;,4e-10
# return (title, pfamid) of the current entry
# >>> ut.blastparse('A0A176VM62_MARPO/110-314 A0A176VM62.1 PF02263.18;GBP;,4e-10')
# 'PF02263'
# used in blasttpfp()
def blastparse(blaststr):
	sarr = blaststr.split(',')
	if len(sarr) < 2:
		cp._err('err:invalid blast str: %s' % blaststr)
	pvalue = float(sarr[1])
	title = sarr[0].split(' ')
	pfam = title[2][0:7]
	if pfam[0:2]!= 'PF':
		cp._err('err:invalid pfam id from: %s' % blaststr)
	return pfam


# return number of tp and fp for one blast PFam output
def blasttpfp(blastoutfile):
	pfamid = blastoutfile[0:7]
	tp = fp = 0
	with open(blastoutfile) as fd:
		for line in fd:
			line = line.strip()
			if len(line) == 0:
				continue
			if blastparse(line) == pfamid:
				tp+=1
			else:
				fp+=1
	if (tp + fp) == 0: # for empty file
		ret = (pfamid, (-1,-1))
	else:
		ret = (pfamid, (tp, fp))
	return ret


# output a column of tp,fp with stub order
def blastcolbystub(arglist):
	if len(arglist) < 3:
		cp._err('Usage: python utils_testcase.py blastcolbystub stubfile blastoutfilelist outfile')

	stubfile = arglist[0]
	ofilelistfile = arglist[1]
	outfile = arglist[2]

	# load stub into a list
	pfamstub = []
	with open(stubfile) as fp:
		for line in fp:
			line = line.strip()
			if len(line)==0:
				continue
			if line[0:2]!='PF':
				cp._err('err:invalid pfam id %s' % line)
			pfamstub.append(line)
	pfamstub.sort() # sort pfamid alphabetically

	# load blast output into a dictionary
	ofilelist = []
	with open(ofilelistfile) as fp:
		for line in fp:
			line = line.strip()
			if len(line) == 0:
				continue
			if line[0:2]!='PF':
				cp._err('err:invalid blast outfile name %s' % line)
			ofilelist.append(line) # PF07582.fa.B62.out

	# for each blast out file calculate tp, fp then store in a dictionary
	# key: pfamid, value: (tp, fp)
	blastout = dict(blasttpfp(file) for file in ofilelist)

	# write outfile with stub order
	with open(outfile, 'w') as fp:
		fp.write('\n'.join(['%s %d %d' % (pfamid, blastout[pfamid][0], blastout[pfamid][1]) for pfamid in pfamstub]))
	cp._info('save to %s' % outfile)


# print tpfp for one blast out file (given expected match)
def blasttpfptuple(arglist):
	if len(arglist) < 2:
		cp._err('Usage: python utils_testcase.py blasttpfptuple outfile match')

	blastoutfile = arglist[0]
	match = arglist[1]
	tp = fp = 0
	with open(blastoutfile) as fd:
		for line in fd:
			line = line.strip()
			if len(line) == 0:
				continue
			if match in line:
				tp+=1
			else:
				fp+=1
	if (tp + fp) == 0: # for empty file
		ret = '-191,-191'
	else:
		ret = '%d %d' % (tp, fp)
	print ret	


# parse FASTA sequence from markx3 align outputs
def parseFasta(lines, i):
	AAset = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X', '-',
			 'a', 'r', 'n', 'd', 'c', 'q', 'e', 'g', 'h', 'i', 'l', 'k', 'm', 'f', 'p', 's', 't', 'w', 'y', 'v', 'b', 'z', 'x']
	fastalines = []
	while lines[i][0] in AAset:
		fastalines.append(lines[i])
		i+=1
	return (''.join(fastalines), i-1)


# parse markx3 format into flat format
def alignparse(title, alignstr):
	lines = filter(None, alignstr.split('\n'))

	i = 0
	out = []
	out.append(title)

	while i<len(lines):
		#print lines[i]
		line = lines[i].strip()
		i+=1

		if 'Program:' in line:
			 tsvline = ' '.join(line.split()) 
			 strArray = tsvline.split(' ')
			 program = strArray[2]
			 out.append(program)
			 #print 'Program: %s' % program
			 continue

		# Matrix: SU
		if 'Matrix:' in line:
			 tsvline = ' '.join(line.split()) 
			 strArray = tsvline.split(' ')
			 matrix = strArray[2]
			 out.append(matrix)
			 #print 'Matrix: %s' % matrix
			 continue

		# Gap_penalty: 10.0
		if 'Gap_penalty:' in line:
			 tsvline = ' '.join(line.split()) 
			 strArray = tsvline.split(' ')
			 gapopen = strArray[2]
			 out.append(gapopen)
			 #print 'Gap_penalty: %s' % gapopen
			 continue

		# Extend_penalty: 0.5			 
		if 'Extend_penalty:' in line:
			 tsvline = ' '.join(line.split()) 
			 strArray = tsvline.split(' ')
			 gapextend = strArray[2]
			 out.append(gapextend)
			 #print 'Extend_penalty: %s' % gapextend
			 continue			

		# Length: 451
		if 'Length:' in line:
			 tsvline = ' '.join(line.split()) 
			 strArray = tsvline.split(' ')
			 length = strArray[2]
			 out.append(length)
			 #print 'length: %s' % length
			 continue

		# Identity:      76/451 (16.9%)	
		if 'Identity:' in line:
			 tsvline = ' '.join(line.split()) 
			 strArray = tsvline.split(' ')
			 ratioStr = strArray[2]
			 ratioArr = ratioStr.split('/')
			 out.append(ratioArr[0]) # number of identity
			 ratio = float(ratioArr[0])/float(ratioArr[1])
			 identity = '%.1f' % (ratio * 100)
			 out.append(identity)
			 #print 'identity: %s' % identity
			 continue

		# Similarity:   125/451 (27.7%)
		if 'Similarity:' in line:
			 tsvline = ' '.join(line.split()) 
			 strArray = tsvline.split(' ')
			 ratioStr = strArray[2]
			 ratioArr = ratioStr.split('/')
			 out.append(ratioArr[0]) # number of similarity
			 ratio = float(ratioArr[0])/float(ratioArr[1])
			 similarity = '%.1f' % (ratio * 100)
			 out.append(similarity)
			 #print 'similarity: %s' % similarity
			 continue

		# Gaps:         186/451 (41.2%)
		if 'Gaps:' in line:
			 tsvline = ' '.join(line.split()) 
			 strArray = tsvline.split(' ')
			 ratioStr = strArray[2]
			 ratioArr = ratioStr.split('/')
			 out.append(ratioArr[0]) # number of gaps
			 ratio = float(ratioArr[0])/float(ratioArr[1])
			 gaps = '%.1f' % (ratio * 100)
			 out.append(gaps)
			 #print 'gaps: %s' % gaps
			 continue

		# Score: 75.5
		if 'Score:' in line:
			 tsvline = ' '.join(line.split()) 
			 strArray = tsvline.split(' ')
			 score = strArray[2]
			 out.append(score)
			 #print 'score: %s' % score
			 continue

		# > ..
		if '>' == line[0]:
			fastaline, index = parseFasta(lines, i) # parse fasta sequence from the ith line
			out.append(('%d' % len(fastaline.replace('-', '')))) # pure sequence length
			out.append(fastaline)
			i = index 
			#print 'seq: %s ... len: %d' % (fastaline[0:5], len(fastaline))
			continue

	return ' '.join(out)


# call needle or water to get the alignment 
# needle <(echo -e ">a\n$1") <(echo -e ">b\n$2") "${@:3}" -filter
def align_exec(seqpair, cmd='needle', matrix='B62', gapopen='10.0', gapextend='0.5'):
	with open(seqpair) as fp:
		line = fp.readline()
		if len(line)<3:
			print 'Error: invalid pairfile : %s\n%s' % (seqpair, line)
			exit(1)
		seqArr = line.split(' ')

	#$ ./align.sh needle "ALIGN" "LINE" B62 8 2
	return sp.Popen(['align.sh', cmd, seqArr[0], seqArr[1], matrix, gapopen, gapextend], stdout=sp.PIPE).communicate()[0]
	#return sp.check_output(['align.sh', cmd, seqArr[0], seqArr[1], matrix, gapopen, gapextend])
 

def printpair(arglist):
	if len(arglist)< 2:
		cp._err('Usage: python proc_testcase.py alignpair pairfile B62')
    
	pairfile = arglist[0]
	matrix = arglist[1]

	ap = palign(alignparse(pairfile, align_exec(pairfile, 'needle', matrix, '10.0', '0.5')))
	ap.dump()


def testpool(arglist):
	if len(arglist)< 5:
		cp.err('Usage: python utils_testcase.py testpool 20p.test.pool needle B62 10.0 0.5')

	poolfile = arglist[0]
	cmd = arglist[1]
	matrix = arglist[2]
	gapopen = arglist[3]
	gapextend = arglist[4]
	outfile = '%s.alnflat' % poolfile

	# loop pool
	identity = 0
	similarity = 0
	gaps = 0

	fout = open(outfile, 'w')
	with open(poolfile) as fp:
		for line in fp:
			name = line.strip()
			ret = align_exec(name, cmd, matrix, gapopen, gapextend) # needle, B62, 10, 0.5
			flat = alignparse(name, ret)
			fout.write('%s\n' % flat)

			ap = palign(flat)
			identity+=ap.nid
			similarity+=ap.nsm
			gaps+=ap.ngp

	fout.close()

	ret =  '%d %d %d %s %s %s %s %s' % (identity, similarity, gaps, poolfile, cmd, matrix, gapopen, gapextend)
	print ret
	return ret


# return number of tp and fp for one blast outfile
# called in tpfpbygap()
def tpfpfromout(blastoutfile):
	# astralS40.00001.a-1-1.fa.b61.6-2.out
	hid = blastoutfile.split('.')[2]
	tp = fp = 0
	with open(blastoutfile) as fd:
		for line in fd:
			# 1aq0A00 3-20-20-80,0.0
			# d2gkma_ a-1-1,3e-89
			line = line.strip()
			if len(line) == 0:
				continue
			s1 = line.split(',')[0]
			outhid = s1.split(' ')[1]
			if outhid == hid:
				tp+=1
			else:
				fp+=1
	return (hid, (tp, fp))


# for coverage vs EPQ test
#				0	   1	 2	 3
# naming code: 	dbname.index.hid.fa
#		 		astralS20.00000.a-1-1.fa
# 				cathS40.00009.3-90-10-10.fa
def tpfpbygap(arglist):
	if len(arglist) < 2:
		cp._err('Usage: python utils_testcase.py tpfpbygap stubfile sm')

	gap = [(9,2), (8,2), (7,2), (6,2), (11,1), (10,1), (9,1)]

	stubfile = arglist[0]
	sm = arglist[1]

	# load homology ID from filename in .stub
	hidset = set()
	namelist = []
	with open(stubfile) as fp:
		for line in fp:
			line = line.strip()
			if len(line)==0:
				continue
			sarr = line.split('.')
			hidset.add(sarr[2])
			namelist.append(line)
	hidlist = list(hidset)
	hidlist.sort()

	# calculate tpfp for different gap configurations
	for g in gap:
		# tpfp dict for all out files with g and sm
		# astralS40.00001.a-1-1.fa.b61.6-2.out
		# name.sm.gap0-gap1.out
		tpfpdict = dict(tpfpfromout('%s.%s.%d-%d.out' % (name, sm, g[0], g[1])) for name in namelist)
		outfile = '%s.%s.%d-%d.tpfp' % (stubfile, sm, g[0], g[1])
		with open(outfile, 'w') as fp:
			fp.write('\n'.join(['%s %d %d' % (hid, tpfpdict[hid][0], tpfpdict[hid][1]) for hid in hidlist]))
			fp.write('\n')
		cp._info('save %s' % outfile)


# generate batch_blast.sh
def batchblastgen(arglist):
	if len(arglist) < 5:
		cp._err('Usage: python utils_testcase2.py batchblastgen dbname stubfile smname format evalue')
	dbname = arglist[0]
	stubfile = arglist[1]
	sm = arglist[2]
	fmt = arglist[3]
	evalue = arglist[4]

	count = 0
	fout = open('batch_blast.sh', 'w')
	querylist = cp.loadlines(stubfile)	
	for seq in querylist:
		for g in cp.gapb62:
			fout.write('blastp -query %s -db $%s -outfmt "%s" -evalue %s -matrix BLOSUM62 -gapopen %d -gapextend %d -out %s.%s.%d-%d.out\n' % (seq, dbname, fmt, evalue, g[0], g[1], seq, sm, g[0], g[1]))
			count+=1
	fout.close()
	cp._info('batch_blast.sh sm: %s len: %d' % (sm, count))

# used in casp target homolog
# re-format .out file into .outflat file
# output: .outflat [target xxxx all --sx-x--x- xxxx B ---dsd---s-s- b62 10-1 1e-4]
# combine .outflat then cut pdb by alignments
# H0953.fa.mc26.sm.9-2.out:
# 	6f45_D,0.0,1,249,MAVQG,1,249,MAVQG
def blastflat_casptarget(arglist):
	if len(arglist) < 2:
		cp._err('Usage: python utils_testcase2.py blastflat_casptarget .tsvfile .outfile')
	tsvfile = arglist[0]
	outfile = arglist[1]

	targetlist = cp.loadlines(tsvfile)
	targetdict = {}
	for target in targetlist:
		sarr = target.split(' ')
		# given target output pdb name
		targetdict[sarr[0]] = sarr[1]

	pdbset = set() # remove all the repeat pdb names xxxx_A,B,C,D
	sarr = outfile.split('.')
	targetname = sarr[0]
	sm = sarr[2]
	gap = sarr[4]
	outlist = cp.loadlines(outfile)

	outfile = outfile + '.outflat'
	fout = open(outfile, 'w')
	for out in outlist:
		sarr = out.split(',')
		uarr = sarr[0].split('_') # 6f45_D
		pdbname = uarr[0]
		if pdbname in pdbset:
			continue
		pdbset.add(pdbname)
		chain = uarr[1]
		evalue = sarr[1]
		qseq = sarr[4]
		sseq = sarr[7]
		outstr = '%s %s all %s %s %s %s %s %s %s\n' % (targetname, targetdict[targetname], qseq, pdbname, chain, sseq, sm, gap, evalue)
		fout.write(outstr)
	fout.close()
	cp._info('save to %s' % outfile)



# combine all .out file into coverage vs EPQ format
def blast2cve(arglist):
	if len(arglist) < 2:
		cp._err('Usage: python utils_testcase2.py blast2cve stubfile sm')

	stubfile = arglist[0]
	sm = arglist[1]

	# load all the fa names
	stublist = []
	with open(stubfile) as fp:
		# naming code: 	dbname.index.hid.fa
		#		 		astralS20.00000.a-1-1.fa
		for faname in fp:
			faname = faname.strip()
			if len(faname)==0:
				continue
			stublist.append(faname)

	if sm not in cp.gapdict:
		cp._info('customized matrix %s, use b62 gaps' % sm)
		gaplist = cp.gapb62
	else:
		gaplist = cp.gapdict[sm]

	#for g in cp.gapb80:
	for g in gaplist:
		outlist = []
		for faname in stublist:
			hid = faname.split('.')[2]
			# load content from .out file
			outname = '%s.%s.%d-%d.out' % (faname, sm, g[0], g[1])
			with open(outname) as fd:
				# d2bkma_ a-1-1,6e-99
				for line in fd:
					line = line.strip()
					if len(line) == 0:
						continue	
					sarr = line.split(',')
					title = sarr[0]
					outhid = title.split(' ')[1]
					tpfp = 1 if hid == outhid else 0
					evalue = float(sarr[1])
					# save tuple	0		1		2	3
					outlist.append((faname,title,evalue,tpfp))
		# sort by evalue
		outlist_sort = sorted(outlist, key=itemgetter(2))
		# save CVE file
		tp=fp=0
		outfile = '%s.%s.%d-%d.cve' % (stubfile, sm, g[0], g[1])
		with open(outfile, 'w') as fout:
			for r in outlist_sort:
				if r[3] == 1:
					tp+=1
				else:
					fp+=1
				outstr = '%s %s %.8f %d %d %d\n' % (r[0], r[1], r[2], r[3], tp, fp)
				fout.write(outstr)
		cp._info('save %s' % outfile)


# output cve points in csv format
# output .tick 
def cvepoint(arglist):
	if len(arglist) < 5:
		cp._info('Usage: python utils_testcase.py cvepoint cvefile nquery epqmax{10} epqnumber outfile')

	cvefile = arglist[0]
	nquery = int(arglist[1])
	mepq = float(arglist[2])
	nepq = int(arglist[3])
	interval = mepq/nepq
	outfile = arglist[4]

	count=0
	cvelist = []
	with open(cvefile) as fp:
		# astralS40.13678.g-52-1.fa d4lgea_ g-52-1 0.00000010 1 71192 24
		for line in fp:
			line = line.strip()
			if len(line)==0:
				continue
			count+=1
			sarr = line.split(' ')
			coverage = float(sarr[5])
			epq = float(sarr[6]) / nquery
			cvelist.append((epq,coverage))
		#cp._info('%d outputs loaded from %s' % (len(cvelist), cvefile))

	# collect coverage on epq tick
	pointer = 0
	cveticks = []
	auc = 0
	for t in np.arange(0, mepq, interval):
		for i in xrange(pointer, len(cvelist)):
			if cvelist[i][0] <= t:
				continue
			else:
				pointer = i
				# tick, coverage, epq
				cveticks.append((t,cvelist[i-1][1], cvelist[i-1][0]))
				auc+=cvelist[i-1][1]
				break

	# write file
	with open(outfile, 'w') as fp:
		fp.write('%s' % ''.join(['%.8f,%.8f,%.8f\n' % (cve[0], cve[1], cve[2]) for cve in cveticks]))
	#cp._info('%d %s %d' % (len(cveticks), outfile, auc))
	print ('%d %s %d' % (len(cveticks), outfile, auc))



##################################################################
# main routine
def main():
	dispatch = {
		'printpair':printpair,
		'testpool':testpool,
		'batchblastgen':batchblastgen,
		'blastflat_casptarget': blastflat_casptarget,
		'blastcolbystub': blastcolbystub,
		'blastcathbystub': blastcathbystub,
		'blasttpfptuple': blasttpfptuple,
		'blast2cve': blast2cve,
		'tpfpbygap': tpfpbygap,
		'cvepoint': cvepoint
	}
	if sys.argv[1] in dispatch:
		dispatch[sys.argv[1]](sys.argv[2:])
	else:
		cp._err('invalid cmd: %s' % sys.argv[1])

if __name__ == '__main__':
	main()
