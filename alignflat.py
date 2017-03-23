import numpy as np
import math
import time
import sys

__all__=['alignflat', 'palign']

'''
$1:  file name
$2:  gap open penalty
$3:  gap extend penalty
$4:  aligned sequence length
$5:  identity number
$6:  identity percentile 
$7:  similarity number
$8:  similarity percentile
$9:  gaps number
$10: gaps percentile
$11: align score
$12: seq A pure length
$13: aligned seq A
$14: seq B pure length
$15: aligned seq B
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

