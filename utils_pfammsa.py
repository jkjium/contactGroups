import sys
import operator
import commp as cp
import numpy as np
import collections


"""
updated utils_msa.py
follow cls -> utils_cls => proc_
"""
class pfammsa(object):
	"""
	data structure of pfam MSA

	"""
	def __init__(self, msafile):
		self.msalist = []
		count = 0
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
	#
	def aafreq(self):
		sumfreq = dict((k,0) for k in cp.msaaa)
		freqlist = [cp.freq(fa[1]) for fa in self.msalist]
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



def aafreq(arglist):
	if len(arglist) < 1 or arglist == False:
		cp._err('Usage: python utils_pfammsa.py aafreq PF00000.txt')

	msafile = arglist[0]

	pfm = pfammsa(msafile)
	output =[(k, float(v)/pfm.seqnum) for k,v in pfm.aafreq().items() if v > 0]
	output.sort(key=operator.itemgetter(1), reverse=True)

	outfile = msafile + '.aafreq'
	with open(outfile, 'w') as fp:
		fp.write(repr(output))
	return cp._info('save to %s' % outfile)


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

	pfm = pfammsa(msafile)
	psubdictall = collections.defaultdict(int)
	for c in scollist:
		psubdict = pfm.pairsubstitution(c[0], c[1])
		for k in psubdict:
			psubdictall[k]+=psubdict[k]

	outfile = '%s.psub' % msafile
	with open(outfile, 'w') as fp:
		for k in psubdictall:
			fp.write('%s %.8f %d %s\n' % (k, float(psubdictall[k])/pfm.msanum, psubdictall[k], cp.quadtype(k)))
	cp._info('save %s' % outfile)

	# return for mp_run reduce
	#return psubdictall


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
	if len(arglist) < 3:
		cp._err('Usage: python utils_pfammsa.py columnselect 1a0p.pdb.A.sgc.cg 1a0p.pdb-PF00589.map PF00589_p90.mip.3.top')

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

	outscolfile = '%s_p90.scol' % (topsdiifile[0:7])
	with open(outscolfile, 'w') as fp:
		fp.write(' '.join([sp for sp in scollist]))

	cp._info('save [%d] selected column(s) to %s' % (len(scollist), outscolfile))


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
		'getsinglemsa': getsinglemsa, # get single MSA gapped / ungapped fa with sequence name or null
		'msareduce': msareduce,
		'pairsubstitution': pairsubstitution,
		'columnselect': columnselect
	}

	if sys.argv[1] not in dispatch:
		print 'invalid cmd: %s' % sys.argv[1]
		return
	else:
		dispatch[sys.argv[1]](sys.argv[2:])

if __name__ == '__main__':
	main()