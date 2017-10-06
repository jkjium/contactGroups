import sys
import operator
import commp as cp
import numpy as np


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

	def aafreq(self):
		sumfreq = dict((k,0) for k in cp.msaaa)
		freqlist = [cp.freq(fa[1]) for fa in self.msalist]
		for dd in freqlist:
			for k in sumfreq:
				sumfreq[k]+=dd[k]
		return sumfreq


def aafreq(arglist):
	if len(arglist) < 1 or arglist == False:
		cp._err('Usage: python utils_pfammsa.py aafreq PF00000.txt')

	msafile = arglist[0]

	pfm = pfammsa(msafile)
	output =[(k, v) for k,v in pfm.aafreq().items() if v > 0]
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



# testing routine
def test():
	pfm = pfammsa('t.pfam.txt')
	pfm.dump()



def main():
	if len(sys.argv)<2:
		print 'Usage: python utils_pfammsa.py cmd pdbfile [args ...]'
		return

	dispatch = {
		'test':test,
		'aafreq': aafreq, # get Amino Acid frequency of a pfam MSA
		'getsinglemsa': getsinglemsa # get single MSA gapped / ungapped fa with sequence name or null
	}

	if sys.argv[1] not in dispatch:
		print 'invalid cmd: %s' % sys.argv[1]
		return
	else:
		dispatch[sys.argv[1]](sys.argv[2:])
	#up = utils_pfammsa('t.pfam.txt')
	#up.dump()
	#up.aafreq()
	#print up.fetchbyname('F7E5S2_MACMU/1044-1120')

if __name__ == '__main__':
	main()