import sys
import operator
import commp as cp
import numpy as np
import collections
import string


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
	def msareduce(self, scoretags, gapcutoff, weightcutoff):
		# scores = { 'aa':[], 'ssp':[] }
		scores = dict((tag, []) for tag in scoretags)
		for s in self.msalist:
			for t in scores:
				scores[t].append([cp.aascore[t][a] for a in s[1]])

		# calculate gap percentage. 0 is reserved for gap
		# generating cp.aascore: //print repr([cp.freq(self.msacol(i))['.']/float(self.msanum) for i in xrange(0, self.msalen)])
		idx_rc = [i for i in xrange(0, self.msalen) if cp.freq(cp.column(scores[scoretags[0]],i))[0]/float(self.msanum) < gapcutoff]

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
		fp.write(repr(output))
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


# generate frequency lookup file for all the column combinations
def freqlookup(arglist):
	if len(arglist) < 3:
		cp._err('Usage: python utils_pfammsa.py freqlookup PF00000_p90.txt.aa.score PF00000_p90.rcol order')

	scorefile = arglist[0]
	colidxfile = arglist[1]
	order = int(arglist[2])
	outfile = '%s.flu.%d' % (scorefile, order)

	data = np.loadtxt(scorefile, delimiter=',')
	colidx = [int(i) for i in np.loadtxt(colidxfile, delimiter=',')]

	with open(outfile, 'w') as fp:
		for s in cp.ncrset(len(colidx), order):
			for (c,lookup) in cp.freqlookup(data[:,s].T):
				#print '%s %s %s %s' % ('-'.join([str(i) for i in s]), ''.join([cp.scoreaa['aa'][c[i]] for i in xrange(order)]), ','.join(['%d' % f for f in c]), ','.join([str(a) for a in lookup]))	
				#print '%s %s %s' % ('-'.join([str(i) for i in s]), ''.join([cp.scoreaa['aa'][c[i]] for i in xrange(order)]), ','.join([str(a) for a in lookup]))	
				fp.write('%s %s %d %s\n' % ('-'.join([str(colidx[i]) for i in s]), ''.join([cp.scoreaa['aa'][c[i]] for i in xrange(order)]), len(lookup), ','.join([str(a) for a in lookup])))
	cp._info('save to %s' % outfile)


def freqlookupscol(arglist):
	if len(arglist) < 3:
		cp._err('Usage: python utils_pfammsa.py freqlookupcol PF00000_p90.txt.aa.score PF00000_p90.rcol PF00000_p90.scol')

	scorefile = arglist[0]
	rcolfile = arglist[1]
	scolfile = arglist[2]
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

	with open(outfile, 'w') as fp:
		for t in scollist:
			s = [colmap[i] for i in t]
			#print repr(t), repr(s)
			for (c,lookup) in cp.freqlookup(data[:,s].T):
				fp.write('%s %s %s %d %s\n' % (scorefile, '-'.join([str(i) for i in t]), ''.join([cp.scoreaa['aa'][c[i]] for i in xrange(order)]), len(lookup), ','.join([str(a) for a in lookup])))
	cp._info('save to %s' % outfile)


# get single MSA gapped / ungapped fa with sequence name or null
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
def msareduce(arglist):
	if len(arglist) < 4:
		cp._err('Usage: $ python utils_pfammsa.py msareduce PF00000.txt aa,ssp 0.7 0.62') 

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
		cp._info('save %s' % outcolfile)

	# output scores
	for t in scoretags:
		outscorefile = '%s.%s.score' % (msafile, t)
		np.savetxt(outscorefile, scores[t], fmt='%d', delimiter=',')
		cp._info('save %s' % outscorefile)


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
		'freqlookup': freqlookup,
		'freqlookupscol': freqlookupscol,
		'getsinglemsa': getsinglemsa, # get single MSA gapped / ungapped fa with sequence name or null
		'msareduce': msareduce,
		'pairsubstitution': pairsubstitution,
		'psicovaln': psicovaln,
		'scoreweight': scoreweight
	}

	if sys.argv[1] not in dispatch:
		print 'invalid cmd: %s' % sys.argv[1]
		return
	else:
		dispatch[sys.argv[1]](sys.argv[2:])

if __name__ == '__main__':
	main()