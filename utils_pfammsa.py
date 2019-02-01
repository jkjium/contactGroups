import sys
import operator
import commp as cp
import numpy as np
import collections
import string
import math


"""
updated utils_msa.py
follow cls -> utils_cls => proc_
"""
class pfammsa(object):
	"""
	data structure of pfam MSA

	"""
	def __init__(self, msafile, opt='full'):
		self.msalist = []
		count = 0

		if opt=='ambaa':
			trans = string.maketrans(''.join(cp.ambaa), ''.join(['.' for i in xrange(len(cp.ambaa))]))
			for head, seq in cp.fasta_iter(msafile):
				self.msalist.append((head, seq.translate(trans)))
		else:
			for head, seq in cp.fasta_iter(msafile):
				#print '%d\n%s\n%s\n' % (count, head, seq)
				self.msalist.append((head, seq))

		self.msalen = len(self.msalist[0][1])
		self.msanum = len(self.msalist)

		if len(self.msalist) == 0:
			cp._err('No MSA found in %s' % msafile)


	def dump(self):
		for s in self.msalist:
			print ">%s\n%s\n" % (s[0], s[1])

	# return ith column AA list
	def msacol(self, i):
		return [s[1][i] for s in self.msalist]

	# return ith column and fa header 
	def msacolfa(self, i):
		return ['>%s\n%s\n' % (s[0], s[1][i]) for s in self.msalist]

	# return a dictionary with dict['A'] = 10
	# cp.freq returns a dictionary for the current sequence
	# input: collist: column index in (int) type, start from 0
	def aafreq(self, collist):
		sumfreq = dict((k,0) for k in cp.msaaa)
		colAAlist = [self.msacol(i) for i in collist]
		freqlist = [cp.freq(c) for c in colAAlist]
		for dd in freqlist:
			for k in sumfreq:
				sumfreq[k]+=dd[k]
		return sumfreq


	# improved msa reduction 
	# converting AA alphabet to aaprop values
	# reduce column by gap percentage
	# assign row weight by hamming distance cluster
	# gapcutoff: gap proportion upper limit
	# weightcutoff is not using. set it -1, external weighting file will be used later
	def msareduce(self, scoretags, gapcutoff, weightcutoff):
		# scores = { 'aa':[], 'ssp':[] }
		scores = dict((tag, []) for tag in scoretags)
		for s in self.msalist:
			for t in scores:
				scores[t].append([cp.aascore[t][a] for a in s[1]])

		# calculate gap percentage. 0 is reserved for gap
		# generating cp.aascore: //print repr([cp.freq(self.msacol(i))['.']/float(self.msanum) for i in xrange(0, self.msalen)])
		idx_rc = [i for i in xrange(0, self.msalen) if cp.freq(cp.column(scores[scoretags[0]],i))[0]/float(self.msanum) < gapcutoff] # allow at most "gapcutoff" % of gap existing in the column

		# output column reduced scores
		scores_rc = {}
		for t in scores:
			scores_rc[t] = np.array(scores[t])[:,idx_rc]

		# weight row
		# add code here

		return scores_rc, idx_rc

	# calculate pair substitution for column i and column j
	def pairsubstitution(self, i, j):
		# frequency of the pairs from two column
		pfreq = cp.freq(['%s%s' % (s[1][i],s[1][j]) for s in self.msalist])
		#print repr(pfreq)
		k = pfreq.keys()
		#print repr(k)
		# count pair substitutions 
		pairsubcount = [('%s%s' % (k[i],k[j]), pfreq[k[i]]*(pfreq[k[j]]-1)/2) if k[i] == k[j] else ('%s%s' % (k[i],k[j]), pfreq[k[i]]*pfreq[k[j]]) for i in xrange(0, len(k)) for j in xrange(i, len(k))]
		#print repr(pairsubcount)

		# combining equivalent quad
		pairsubdict = collections.defaultdict(int)
		for quad, count in pairsubcount:
			pairsubdict[cp.quad_permu(list(quad))]+=count

		return pairsubdict


# calculate AA frequency by given column indices
def aafreq(arglist):
	if len(arglist) < 1 or arglist == False:
		cp._err('Usage: python utils_pfammsa.py aafreq PF00000.txt col=all|1,2,3')

	msafile = arglist[0]
	pfm = pfammsa(msafile)
	if arglist[1] == 'all':
		collist = [i for i in xrange(pfm.msalen)]
	else:
		collist = [int(i) for i in arglist[1].split(',')]
	print repr(collist)

	freqdict = pfm.aafreq(collist)
	freqsum = sum(freqdict.values())
	output =[(k, v, float(v)/freqsum) for k,v in freqdict.items() if v > 0]
	output.sort(key=operator.itemgetter(1), reverse=True)

	outfile = msafile + '.aafreq'
	with open(outfile, 'w') as fp:
		fp.write('%s\n' % (','.join(['%s %d %.8f' % (k, v, nv) for (k,v,nv) in output])))
	return cp._info('save to %s' % outfile)


# calculate AA frequency by given selected column file
def aafreqscol(arglist):
	if len(arglist) < 2 or arglist == False:
		cp._err('Usage: python utils_pfammsa.py aafreq PF00000.txt PF00000_p90.scol')

	msafile = arglist[0]
	scolfile = arglist[1]

	# load scols
	scolset = set()
	with open(scolfile) as fp:
		# 26-82 34-76
		for line in fp:
			line = line.strip()
			if len(line) == 0:
				continue
			for colpair in line.split(' '):
				colarr = colpair.split('-')
				scolset.add(int(colarr[0]))
				scolset.add(int(colarr[1]))
	#print repr(scolset)

	pfm = pfammsa(msafile)
	freqdict = pfm.aafreq(scolset)
	freqsum = sum(freqdict.values())
	output =[(k, v, float(v)/freqsum) for k,v in freqdict.items() if v > 0]
	output.sort(key=operator.itemgetter(1), reverse=True)

	outfile = msafile + '.aafreqscol' 
	with open(outfile, 'w') as fp:
		fp.write('%s\n' % (','.join(['%s %d %.8f' % (k, v, nv) for (k,v,nv) in output])))
	return cp._info('save to %s' % outfile)



# select key columns by residue contact & top sdii
# input: 1a0p.pdb.A.sgc.cg, 1a0p.pdb-PF00589.map, PF00589.mip.3.top
# $ head PF00589_p90.mip.3.top
#	271-274 0.25128949
# $ head 1a0p.pdb.A.sgc.cg
#	3 Q 7 R
# $ tail 1a0p.pdb-PF00589.map
#	280 T 1318 I
# output: PF00589_p90.scol
def columnselect(arglist):
	if len(arglist) < 4:
		cp._err('Usage: python utils_pfammsa.py columnselect 1a0p.pdb.A.sgc.cg 1a0p.pdb-PF00589.map PF00589_p90.mip.3.top PF00589_p90_tip.scol')

	cgfile = arglist[0]
	with open(cgfile) as fp:
		# [['3', 'Q', '7', 'R'], ['4', 'D', '7', 'R']...]
		cglist = [line.strip().split(' ') for line in fp]
	#print repr(cglist)

	mapfile = arglist[1]
	with open(mapfile) as fp:
		maplist = [line.strip().split(' ') for line in fp]
	#print repr(maplist)

	topsdiifile = arglist[2]
	with open(topsdiifile) as fp:
		topsdiiset = set([line.strip().split(' ')[0] for line in fp])
	#print repr(topsdiiset)

	# convert maplist to resid -> msa dictionary
	mapdict = dict((m[0], m[2]) for m in maplist)
	#print repr(mapdict)
	
	scollist = []
	# choose pair that in contact and top sdii	
	for c in cglist:
		# c: ['113', 'L', '170', 'I']
		# contact resi not in mapdict
		if (c[0] not in mapdict) or (c[2] not in mapdict):
			continue
		p1 = int(mapdict[c[0]])
		p2 = int(mapdict[c[2]])
		key = '%d-%d' % (p1, p2) if p1 <= p2 else '%d-%d' % (p2, p1)
		if key in topsdiiset:
			scollist.append(key)

	if len(scollist) == 0:
		cp._info('err:no significant columnt selected for %s' % topsdiifile[0:7])
		return

	outscolfile = arglist[3]
	with open(outscolfile, 'w') as fp:
		fp.write(' '.join([sp for sp in scollist]))

	cp._info('save [%d] selected column(s) to %s' % (len(scollist), outscolfile))


# for homolog testcase
# save the first entry to a fasta file
def seqfa(arglist):
	if len(arglist) < 2:
		cp._err('Usage: python utils_pfammsa.py msafile seqfa outfile')

	msafile = arglist[0]
	outfile = arglist[1]

	with open(outfile, 'w') as fp:
		for head, msa in cp.fasta_iter(msafile):
			seq = msa.translate(None, ''.join(cp.gaps))
			fp.write('>%s\n%s' % (head, seq))
			cp._info('save to %s' % outfile)
			break # only get the first one



# calculate mean entropy for score file
def scoreentropy(arglist):
	if len(arglist) < 1:
		cp._err('Usage: python utils_pfammsa.py scoreentropy scorefile')

	scorefile = arglist[0]
	score = np.loadtxt(scorefile, delimiter=',')
	hlist = [cp.entropy([score[:,i]]) for i in xrange(0, score.shape[1])]
	print '%s %.8f' % (scorefile, sum(hlist)/score.shape[1])


# generate frequency lookup file for all the column combinations
def freqlookup(arglist):
	if len(arglist) < 3:
		cp._err('Usage: python utils_pfammsa.py freqlookup PF00000_p90.txt.aa.score PF00000_p90.rcol order')

	scorefile = arglist[0]
	colidxfile = arglist[1]
	order = int(arglist[2])
	outfile = '%s.%d.flu' % (scorefile, order)

	data = np.loadtxt(scorefile, delimiter=',')
	colidx = [int(i) for i in np.loadtxt(colidxfile, delimiter=',')]

	with open(outfile, 'w') as fp:
		for s in cp.ncrset(len(colidx), order):
			for (c,lookup) in cp.freqlookup(data[:,s].T):
				#print '%s %s %s %s' % ('-'.join([str(i) for i in s]), ''.join([cp.scoreaa['aa'][c[i]] for i in xrange(order)]), ','.join(['%d' % f for f in c]), ','.join([str(a) for a in lookup]))	
				#print '%s %s %s' % ('-'.join([str(i) for i in s]), ''.join([cp.scoreaa['aa'][c[i]] for i in xrange(order)]), ','.join([str(a) for a in lookup]))	
				fp.write('%s %s %.4f %s\n' % ('-'.join([str(colidx[i]) for i in s]), ''.join([cp.scoreaa['aa'][c[i]] for i in xrange(order)]), float(len(lookup))/len(data), ','.join([str(a) for a in lookup])))
	cp._info('save to %s' % outfile)


# calculate freq lookup table for each selected column per pfam
# scol is filtered by map in column selection function
def freqlookupscol(arglist):
	if len(arglist) < 4:
		cp._err('Usage: python utils_pfammsa.py freqlookupscol PF00000_p90.txt.aa.score PF00000_p90.rcol PF00000_p90.scol 1234.pdb-PF00000_p90.map')

	scorefile = arglist[0]
	rcolfile = arglist[1]
	scolfile = arglist[2]
	rmapfile = arglist[3]
	outfile = '%s.scolflu' % scorefile

	# load score
	data = np.loadtxt(scorefile, delimiter=',')
	rcol = np.loadtxt(rcolfile, delimiter=',')
	colmap = dict((int(rcol[i]), i) for i in xrange(len(rcol)))

	# load column tuple
	# 274-474 359-366 386-400
	scollist = []
	with open(scolfile) as fp:
		line = fp.readline().strip()
		if len(line)!=0:
			pairstr = line.split()
			for p in pairstr:
				sarr = p.split('-')
				order = len(sarr)
				scollist.append([int(c) for c in sarr])

	# load map 
	# (resi  msai)  (279 Y 1317 Y)
	msai2resi = {}
	with open(rmapfile) as fp:
		for line in fp:
			line = line.strip()
			if len(line)!=0:
				sarr = line.split(' ')
				msai2resi[int(sarr[2])] = int(sarr[0])
	#print repr(msai2resi)

	# write lookup
	with open(outfile, 'w') as fp:
		for t in scollist:
			s = [colmap[i] for i in t]
			#print repr(t), repr(s)
			for (c,lookup) in cp.freqlookup(data[:,s].T):
				fp.write('%s %s %s %s %.4f %s\n' % (scorefile, '-'.join([str(i) for i in t]), '-'.join([str(msai2resi[i]) for i in t]),
												''.join([cp.scoreaa['aa'][c[i]] for i in xrange(order)]), 
												float(len(lookup))/len(data), 
												','.join([str(a) for a in lookup])))

	cp._info('save to %s' % outfile)

# get a single column from MSA
def getcolumn(arglist):
	if len(arglist) < 3:
		cp._err('Usage: python utils_pfammsa.py getcolumn pfamID column_idx{start from 0} outfile')

	msafile = arglist[0]
	colindex = int(arglist[1])
	outfile = arglist[2]

	pfm = pfammsa(msafile)
	if colindex >= pfm.msalen:
		cp._err('colindex : %d is larger than MSA len %d' % (colindex, pfm.msalen))
	colstrlist = pfm.msacolfa(colindex)
	with open(outfile, 'w') as fout:
		fout.write(''.join(colstrlist))
	cp._info('column data save to %s' % outfile)


# get single MSA gapped / ungapped fa with sequence name or null
'''
output: save [null] raw seq : PF00008_p90_seq.fa , MSA seq : PF00008_p90_MSA.fa
PF00008_p90_seq.fa: ungapped fa, for pfamscan
PF00008_p90_MSA.fa: gapped fa, for resimap
'''
def getsinglemsa(arglist):
	if len(arglist) < 1:
		cp._err('Usage: python utils_pfammsa.py getsinglemsa PF00000.txt head PF00000')

	msafile = arglist[0]
	head = arglist[1]
	outprefix = arglist[2]

	pfm = pfammsa(msafile)
	# with head
	if head!='null':
		for s in pfm.msalist:
			if head in s:
				msa = s
	# no head, get the first MSA
	else:
		msa = pfm.msalist[0]

	outseqfile = '%s_seq.fa' % outprefix
	with open(outseqfile, 'w') as fp:
		fp.write('>%s\n%s' % (msa[0],msa[1].translate(None, ''.join(cp.gaps))))

	outmsafile = '%s_MSA.fa' % outprefix
	with open(outmsafile, 'w') as fp:
		fp.write('>%s\n%s' % (msa[0],msa[1]))

	cp._info('save [%s] raw seq : %s , MSA seq : %s' % (head, outseqfile, outmsafile))



# improved version of utils_msa.MSAReduction()
# input : PF0000.txt, scoretags, gapcutoff(max gap percentage), weightcutoff
# output: PF0000.txt.{stag.score, col, row}
# python utils_pfammsa.py msareduce t.txt aa 0.2 -1 // save columns have < 20% of gaps
def msareduce(arglist):
	if len(arglist) < 4:
		cp._err('Usage: $ python utils_pfammsa.py msareduce PF00000.txt aa,ssp 0.7 -1') 

	msafile = arglist[0]
	scoretags = arglist[1].split(',')
	gapcutoff = float(arglist[2])
	weightcutoff = float(arglist[3])

	pfm = pfammsa(msafile)
	#pfm.dump()
	scores, idx_rc = pfm.msareduce(scoretags, gapcutoff, weightcutoff)

	# output column
	outcolfile = msafile + '.rcol'
	with open(outcolfile, 'w') as fp:
		fp.write(','.join([str(i) for i in idx_rc]))
		cp._info('save %d/%d columns in %s' % (len(idx_rc), pfm.msalen, outcolfile))

	# output scores
	for t in scoretags:
		outscorefile = '%s.%s.score' % (msafile, t)
		np.savetxt(outscorefile, scores[t], fmt='%d', delimiter=',')
		cp._info('save %d rows in %s' % (len(scores[t]), outscorefile))


# calculate the pair substitution for one PFam MSA
# input: PF00000.txt, PF00000.txt.selected.col
# output: tsv file with quad subsitution count
def pairsubstitution(arglist):
	if len(arglist) < 2:
		cp._err('Usage: python utils_pfammsa.py pairsubstitution PF00000.txt PF00000.txt.scol')

	msafile = arglist[0]
	scolfile = arglist[1]

	# load selected columns 
	# 1,2 5,8 3,10
	with open(scolfile) as fp:
		line = fp.readline().strip()
		if len(line) == 0:
			cp._info('no column loaded for %s, %s' % (msafile, scolfile))
			return
		else:
			scollist = [[int(c) for c in p.split('-')] for p in line.split(' ')]

	# get scol set from pairs list
	scolset = set()
	for p in scollist:
		scolset.add(p[0])
		scolset.add(p[1])

	#pfm = pfammsa(msafile)
	pfm = pfammsa(msafile, 'ambaa')
	fq = pfm.aafreq(scolset)
	psubdictall = collections.defaultdict(int)
	# counting freq
	for c in scollist:
		psubdict = pfm.pairsubstitution(c[0], c[1])
		for k in psubdict:
			psubdictall[k]+=psubdict[k]

	# normalize by AA freq
	# k[0] k[1] k[2] k[3] must in freqdict
	norm_psubdictall = cp.rank01(dict((k, float(psubdictall[k])/(fq[k[0]]*fq[k[1]]*fq[k[2]]*fq[k[3]])) for k in psubdictall))

	outfile = '%s.psub' % msafile
	with open(outfile, 'w') as fp:
		for k in psubdictall:
			fp.write('%s %s %d %.8f %.8f\n' % (k, cp.quadtype(k), psubdictall[k], float(psubdictall[k])/pfm.msanum, norm_psubdictall[k]))
	cp._info('save %s' % outfile)

	# return for mp_run reduce
	#return psubdictall


def psicovaln(arglist):
	if len(arglist) < 1:
		cp._err('Usage: python utils_pfammsa.py psicovaln PF00000_full.txt')
	msafile = arglist[0]
	outfile = msafile + '.aln'
	pfm = pfammsa(msafile,'ambaa')
	with open(outfile, 'w') as fp:
		for s in pfm.msalist:
			fp.write('%s\n' % s[1].replace('.', '-'))
	cp._info('save to %s' % outfile)


# print scol -> resi set
def scol2resi(arglist):
	if len(arglist) < 2:
		cp._err('Usage: python utils_pfammsa.py scol2resi mapfile scolfile')

	mapfile = arglist[0]
	scolfile =arglist[1]

	resimap = {}
	with open(mapfile) as fp:
		for line in fp:
			line = line.strip()
			if len(line)==0:
				continue
			# 7 W 52 I
			sarr = line.split(' ')
			resimap[sarr[2]] = sarr[0]

	scolset = set()
	with open(scolfile) as fp:
		line = fp.readline()
		line = line.strip()
		if len(line) != 0:
			for p in line.split(' '):
				sarr = p.split('-')
				scolset.add(sarr[0])
				scolset.add(sarr[1])

	if len(scolset) == 0:
		print '%s -1' % (mapfile)
	else:
		resi = [int(resimap[c]) for c in scolset]
		resi.sort()
		print '%s %s' % (mapfile, ' '.join([str(r) for r in resi]))




# calculate and save sequence weight with similarity cutoff
def scoreweight(arglist):
	if len(sys.argv) < 2:
		cp._err('Usage: python utils_pfammsa.py scoreweight PF00000.score similarity_value')

	datafile = arglist[0]
	svalue = float(arglist[1])
	outfile = '%s.%2d.w' % (datafile, svalue*100)

	score = np.loadtxt(datafile, delimiter=',')
	w = cp.hamming_weight(score, 1-svalue)
	np.savetxt(outfile, w)
	cp._info('save weight to %s' % outfile)


# split fasta file into separate .fa file
# filename: prefix.00001.fa
def splitfa(arglist):
	if len(arglist)<3:
		cp._err('Usage: python utils_pfammsa.py splitfa fastfile outprefix padding')

	fastafile = arglist[0]
	outprefix = arglist[1]
	padding = int(arglist[2])

	count=0
	for head, seq in cp.fasta_iter(fastafile):
		outfafile = '%s.%s.fa' % (outprefix, str(count).zfill(padding))
		with open(outfafile, 'w') as fp:
			fp.write('>%s\n%s' % (head, seq))
		count+=1
	cp._info('%s : save %d .fa files' % (fastafile, count))


# for coverage vs EPQ homology testing
# split fasta file into separate .fa file
# filename: prefix.00001.a-1-1.fa
def splithidfa(arglist):
	if len(arglist)<3:
		cp._err('Usage: python utils_pfammsa.py splithidfa fastfile outprefix padding')

	fastafile = arglist[0]
	outprefix = arglist[1]
	padding = int(arglist[2])

	count=0
	for head, seq in cp.fasta_iter(fastafile):
		hid = head.split(' ')[1]
		outfafile = '%s.%s.%s.fa' % (outprefix, str(count).zfill(padding), hid)
		with open(outfafile, 'w') as fp:
			fp.write('>%s\n%s' % (head, seq))
		count+=1
	cp._info('%s : save %d .fa files' % (fastafile, count))



# accumulate weighted background frequcney for each AA 
def wfreq(arglist):
	if len(arglist) < 4:
		cp._err('Usage: python utils_pfammsa.py wfreq PF01012_p90.txt.score.flu.1 PF00107_p90.txt.62.weight PF13507_p90_tip.scol PF01012.w70.wfreq')

	flufile = arglist[0]
	wfile = arglist[1]
	scolfile = arglist[2]
	outfile = arglist[3]

	# load w
	if wfile != 'na':
		w = np.loadtxt(wfile)

	# load selected col
	scolset = set()
	with open(scolfile) as fp:
		line = fp.readline().strip()
		if len(line.strip())==0:
			cp._info('%s no scol pairs avaiable' % scolfile)
			return
		for p in line.split(' '):
			c = p.split('-')
			scolset.add(c[0])
			scolset.add(c[1])
	#print repr(scolset)		

	wfreqdict = collections.defaultdict(float)	 # for denominator from all column
	wsfreqdict = collections.defaultdict(float)	 # for denominator from scol
	wscoldict =	collections.defaultdict(list) 
	# calculate sinlgle wfreq
	# 206 C .5,111,133,158,210,241
	with open(flufile) as fp:
		for line in fp:
			line = line.strip()
			if len(line) == 0:
				continue
			sarr = line.split(' ')
			if wfile == 'na':
				wfreq = len(sarr[3].split(','))
			else:
				wfreq = sum(w[[int(i) for i in sarr[3].split(',')]])
			wfreqdict[sarr[1]]+= wfreq
			#print '%s, %s, %s, %.4f' % (sarr[1], [int(i) for i in sarr[3].split(',')], repr(w[[int(i) for i in sarr[3].split(',')]]), wfreqdict[sarr[1]])
			# save freq of scol for substitution frequency calculation
			if sarr[0] in scolset:
				wscoldict[sarr[0]].append((sarr[1], wfreq))
				wsfreqdict[sarr[1]]+= wfreq
	#print repr(wscoldict)

	# calculate substitution wfreq
	wsmdict = collections.defaultdict(float)
	for c in wscoldict: # {'205': [('C', 0.7), ('D', 0.5)], '207': [('C', 1.0)]}
		wfl = wscoldict[c]
		for i in xrange(len(wfl)):
			A = wfl[i][0]
			if A in cp.abaa:
				continue
			for j in xrange(i+1, len(wfl)):
				B = wfl[j][0]
				if B in cp.abaa:
					continue
				k = A+B if A < B else B+A
				wsmdict[k]+=wfl[i][1]*wfl[j][1]
			wsmdict[A+A]+=wfl[i][1]*wfl[i][1]/2

	# output
	with open(outfile, 'w') as fp:
		# save single wfreq
		for c in wfreqdict:
			fp.write('wf %s %.8f\n' % (c, wfreqdict[c]))
		# save single wfreq from scol
		for c in wsfreqdict:
			fp.write('sf %s %.8f\n' % (c, wsfreqdict[c]))
		# save substitution wfreq
		for k in wsmdict:
			fp.write('sd %s %.8f\n' %(k, wsmdict[k]))
	cp._info('save wfreq to %s' % outfile)



# accumulate weighted background frequcney for each AA 
# from sscol: 23 32 151
def wfreqs(arglist):
	if len(arglist) < 4:
		cp._err('Usage: python utils_pfammsa.py wfreqs PF01012_p90.txt.score.flu.1 PF00107_p90.txt.62.weight PF13507_p90_tip.sscol PF01012.w70.wfreq')

	flufile = arglist[0]
	wfile = arglist[1]
	scolfile = arglist[2]
	outfile = arglist[3]

	# load w
	if wfile != 'na':
		w = np.loadtxt(wfile)

	# load selected col
	scolset = set()
	with open(scolfile) as fp:
		line = fp.readline().strip()
		if len(line.strip())==0:
			cp._info('%s no scol pairs avaiable' % scolfile)
			return
		scolset = set([p for p in line.split(',')])
	#print repr(scolset)		
	
	wfreqdict = collections.defaultdict(float)	 # for denominator from all column
	wsfreqdict = collections.defaultdict(float)	 # for denominator from scol
	wscoldict =	collections.defaultdict(list) 
	# calculate sinlgle wfreq
	# 206 C .5,111,133,158,210,241
	with open(flufile) as fp:
		for line in fp:
			line = line.strip()
			if len(line) == 0:
				continue
			sarr = line.split(' ')
			if wfile == 'na':
				wfreq = len(sarr[3].split(','))
			else:
				wfreq = sum(w[[int(i) for i in sarr[3].split(',')]])
			wfreqdict[sarr[1]]+= wfreq
			#print '%s, %s, %s, %.4f' % (sarr[1], [int(i) for i in sarr[3].split(',')], repr(w[[int(i) for i in sarr[3].split(',')]]), wfreqdict[sarr[1]])
			# save freq of scol for substitution frequency calculation
			if sarr[0] in scolset:
				wscoldict[sarr[0]].append((sarr[1], wfreq))
				wsfreqdict[sarr[1]]+= wfreq
	#print repr(wscoldict)

	# calculate substitution wfreq
	wsmdict = collections.defaultdict(float)
	for c in wscoldict: # {'205': [('C', 0.7), ('D', 0.5)], '207': [('C', 1.0)]}
		wfl = wscoldict[c]
		for i in xrange(len(wfl)):
			A = wfl[i][0]
			if A in cp.abaa:
				continue
			for j in xrange(i+1, len(wfl)):
				B = wfl[j][0]
				if B in cp.abaa:
					continue
				k = A+B if A < B else B+A
				wsmdict[k]+=wfl[i][1]*wfl[j][1]
			wsmdict[A+A]+=wfl[i][1]*wfl[i][1]/2

	# output
	with open(outfile, 'w') as fp:
		# save single wfreq
		for c in wfreqdict:
			fp.write('wf %s %.8f\n' % (c, wfreqdict[c]))
		# save single wfreq from scol
		for c in wsfreqdict:
			fp.write('sf %s %.8f\n' % (c, wsfreqdict[c]))
		# save substitution wfreq
		for k in wsmdict:
			fp.write('sd %s %.8f\n' %(k, wsmdict[k]))
	cp._info('save wfreq to %s' % outfile)


# combine single frequency and substitution frequency into sm
def wfreq2sm(arglist):
		if len(arglist) < 3:
			cp._err('Usage: python utils_pfammsa.py wfreq2sm combine.wfreq wf|sf outfile')

		wfreqfile = arglist[0]
		opt = arglist[1]
		outprefix = arglist[2]

		qij = collections.defaultdict(float)
		eij = collections.defaultdict(float)
		with open(wfreqfile) as fp:
				for line in fp:
						line = line.strip()
						if len(line)==0:
								continue
						sarr = line.split(' ')
						t = sarr[0]
						k = sarr[1]
						f = float(sarr[2])
						# sf M 3511.29
						if t == opt:
						#if len(k) == 1:
								eij[k]+=f
						# sd QR 95488.206
						elif t == 'sd':
						#elif len(k)== 2:
								qij[k]+=f

		# single freq
		total_e = sum(eij.values())
		for k in eij:
				eij[k]=eij[k]/total_e
		# substitution freq
		total_q = sum(qij.values())
		for k in qij:
				qij[k]=qij[k]/total_q

		#print 'len(qij): %d' % len(qij)
		# calculate sm
		sm = collections.defaultdict(int)
		for k in qij:
			A = k[0]
			B = k[1]
			if A==B:
				sm[A+B] = int(round(2*math.log(qij[A+B]/(eij[A]*eij[B]),2)))
			else:
				sm[A+B] = int(round(2*math.log(qij[A+B]/(2*eij[A]*eij[B]),2)))
			sm[B+A] = sm[A+B]
		#print min(sm.values()), max(sm.values())

		# save raw sm
		'''
		with open(outprefix+'.raw.sm', 'w') as fp:
			for A in cp.aat01:
				fp.write('%s\n' % ' '.join([str(sm[A+B]).rjust(2) for B in cp.aat01]))
		cp._info('save raw sm to %s' % outprefix)

		# output readable sm
		npstd = np.array([[sm[A+B] for B in cp.aat01] for A in cp.aat01])
		stdfile = outprefix + '.std.sm'
		with open(stdfile, 'w') as fp:
			fp.write(cp.smstr(npstd, cp.aat01))
		cp._info('save std sm to %s' % stdfile)
		'''
		
		# output emboss sm
		embossfile = outprefix + '.emboss.sm'
		npemboss = np.array([[sm[A+B] for B in cp.smaa2] for A in cp.smaa2])
		cp.b62edge[:npemboss.shape[0], :npemboss.shape[1]] = npemboss
		with open(embossfile, 'w') as fp:
			fp.write(cp.smstr(cp.b62edge, cp.smaa1))
		cp._info('save emboss sm to %s' % embossfile)


# testing routine
def test(arglist):
	# test columnselect

	'''
	>A2SSP0_METLZ/1-86
	MH..E..F

	>A3CU99_METMJ/1-86
	MH..E..D

	>A5UJB7_METS3/1-88
	MY..A.AC
	'''
	'''
	pfm = pfammsa('PF00000.txt')
	pfm.dump()
	ret = pfm.pairsubstitution(0,1)
	print repr(ret)
	'''
	# test pfm.msareduce()
	'''
	pfm = pfammsa('PF00000.txt')
	pfm.dump()
	#print pfm.msacol(1)
	scores, idx_rc = pfm.msareduce(['aa', 'ssp'],0.7,0.0)
	print repr(idx_rc)
	for sc in scores:
		print sc
		print repr(scores[sc])
	'''


def main():
	if len(sys.argv)<2:
		print 'Usage: python utils_pfammsa.py cmd pdbfile [args ...]'
		return

	dispatch = {
		'test':test,
		'aafreq': aafreq, # get Amino Acid frequency of a pfam MSA
		'aafreqscol': aafreqscol,
		'columnselect': columnselect,
		'seqfa': seqfa,
		'scoreentropy': scoreentropy,
		'freqlookup': freqlookup,
		'freqlookupscol': freqlookupscol,
		'getcolumn': getcolumn, # get a single column from MSA
		'getsinglemsa': getsinglemsa, # get single MSA gapped / ungapped fa with sequence name or null
		'msareduce': msareduce,
		'pairsubstitution': pairsubstitution,
		'psicovaln': psicovaln,
		'scol2resi': scol2resi,
		'scoreweight': scoreweight,
		'splitfa': splitfa,
		'splithidfa': splithidfa,
		'wfreq': wfreq,
		'wfreqs': wfreqs,
		'wfreq2sm': wfreq2sm
	}

	if sys.argv[1] not in dispatch:
		print 'invalid cmd: %s' % sys.argv[1]
		return
	else:
		dispatch[sys.argv[1]](sys.argv[2:])

if __name__ == '__main__':
	main()
