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

		self.scoreBinary = {
							'-': 0,'.': 0,'A': 1,'C': 1,'D': 1,'E': 1,'F': 1,'G': 1,'H': 1,'I': 1,'K': 1,
							'L': 1,'M': 1,'N': 1,'P': 1,'Q': 1,'R': 1,'S': 1,'T': 1,'V': 1,'W': 1,'Y': 1
						}

		self.seqlen = len(self.msaArray[0][1])
		self.seqNum = len(self.msaArray)


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


	# generate concised msa scoreboard
	def msaboard(self, cutoff):
		scoreboard = []
		addboard = np.zeros(self.seqlen)
		for s in self.msaArray:
			scoreboard.append([self.scoreValue[a.upper()] for a in s[1]])
			addboard+=np.array([self.scoreBinary[a.upper()] for a in s[1]]) # calculate gap proportion

		# get conserved columns
		indices = [i for i in xrange(0, self.seqlen) if addboard[i]/self.seqNum  > cutoff]
		#print np.array(scoreboard)
		return (np.array(scoreboard)[:,indices], indices)


	# write scoreboard
	def writeScoreboard(self, outfile, cutoff):
		fout = open(outfile, 'w')
		score, varlist = self.msaboard(cutoff)
		fout.write('%s\n' % repr(varlist))
		for line in score:
			fout.write('%s\n' % ','.join([str(i) for i in line]))
		fout.close()

