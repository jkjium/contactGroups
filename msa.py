import sys
import numpy as np
from itertools import groupby

class msa(object):
	def __init__(self, msafile):
		self.msaArray=[]
		for s in self.fasta_iter(msafile):
			self.msaArray.append(s)

		self.scoreValue = {
							'-': 0,'.': 0,'A': 1,'C': 2,'D': 3,'E': 4,'F': 5,'G': 6,'H': 7,'I': 8,'K': 9,
							'L': 10,'M': 11,'N': 12,'P': 13,'Q': 14,'R': 15,'S': 16,'T': 17,'V': 18,'W': 19,'Y': 20
						}

	# given a fasta file yield tuples of header, sequence
	def fasta_iter(self, fastafile):
		fh = open(fastafile)
		faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))	
		for header in faiter:
			header = header.next()[1:].strip()
			seq = "".join(s.strip() for s in faiter.next())
			yield header, seq


	def dump(self):
		for s in self.msaArray:
			print s[0]
			print s[1]
			print
		print '%d sequence in total' % len(self.msaArray)


	# get corresponding map from sequence postion to msa position
	def getPosMap(self, p):
		pdbseq = p.seq
		pdbheader = p.pdb

		msaheader = ''
		msaseq = ''
		for s in self.msaArray:
			msaheader = s[0]
			if pdbheader in msaheader:
				msaseq = s[1]
				break
		if msaseq == '':
			print 'getPosMap()::error:header %s not found' % pdbheader
			return

		idxMap = {}
		token = -1
		for i in xrange(0,len(pdbseq)):
			for j in xrange(token+1, len(msaseq)):
				if pdbseq[i] == msaseq[j]:
					idxMap[i] = j
					token = j
					break

		if len(idxMap) != len(pdbseq):
			print 'getPosMap()::error:length not match'
			return

		return idxMap


	# generate msa scoreboard
	def msaboard(self):
		scoreboard = []
		for s in self.msaArray:
			scoreboard.append([self.scoreValue[a.upper()] for a in s[1]])

		return scoreboard

	# write scoreboard
	def writeScoreboard(self, outfile):
		fout = open(outfile, 'w')
		score = self.msaboard()
		for line in score:
			fout.write('%s\n' % ','.join([str(i) for i in line]))
		fout.close()

