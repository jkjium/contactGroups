import sys
import numpy as np
from itertools import groupby

class msa(object):
	def __init__(self, msafile):
		self.target = ('', '')
		self.msaArray=[]
		for s in self.fasta_iter(msafile):
			self.msaArray.append(s)
		'''	
			if target in s[0]:
				self.target = s
		if self.target[1] == '':
			print 'Target: %s not found in MSA.' % target
			sys.exit(-1)
		'''
		self.scoreValue = {
							'X':0,'-': 0,'.': 0,'A': 1,'C': 2,'D': 3,'E': 4,'F': 5,'G': 6,'H': 7,'I': 8,'K': 9,
							'L': 10,'M': 11,'N': 12,'P': 13,'Q': 14,'R': 15,'S': 16,'T': 17,'V': 18,'W': 19,'Y': 20, 'B': 3
						}

		self.scoreBinary = {
							'X':0,'-': 0,'.': 0,'A': 1,'C': 1,'D': 1,'E': 1,'F': 1,'G': 1,'H': 1,'I': 1,'K': 1,
							'L': 1,'M': 1,'N': 1,'P': 1,'Q': 1,'R': 1,'S': 1,'T': 1,'V': 1,'W': 1,'Y': 1, 'B': 1
						}

		self.seqlen = len(self.msaArray[0][1])
		self.seqNum = len(self.msaArray)

		# sequence weight in msa
		self.weight = np.ones(self.seqNum)

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


	# given pdb sequence and target in msa
	# return pdb resi -> msa position (with gaps) map
	# for pdb sequence contained in the msa
	def getResiTargetMap(self, p, target):
		midict = {} # midict[resn+seq_pos] = msa_pos 
		rtdict = {}
		targetSeqArr = []

		pdbseq = p.seq
		self.setTarget(target)

		targetmsa = self.target[1].upper()

		index = 0
		for i in xrange(0, len(targetmsa)):
			a = targetmsa[i]
			if a not in ['-', '.']:
				targetSeqArr.append(a)
				key = '%s%d' % (a, index) # resn with seq index
				midict[key] = i # i is the msa index
				index+=1

		targetseq = ''.join(targetSeqArr)

		print 'pdb seq:\n%s' % (pdbseq)
		print 'target seq:\n%s' % (targetseq)
		preffixp = pdbseq.find(targetseq)
		if preffixp < 0:
			print 'Cannot match target sequence with pdb sequence'
			print 'pdb seq:\n%s' % (pdbseq)
			print 'target seq:\n%s' % (targetseq)

		# p.resDict : (index in sequence, ResName)
		# 'B529': (132, 'V')
		for k in p.resDict:
			(pos, resn) = p.resDict[k]
			resn = resn.upper()
			mikey = '%s%d' % (resn, pos - preffixp)
			if mikey in midict:
 				msapos = midict[mikey]
 			else:
 				continue # skip the key like (k-2)

			if targetmsa[msapos] != resn:
				print 'msa:getResiTargetMap():: mismatch from pdb to msa'
				print 'msa: [%s], pdb: [%s]' % (self.target[1][msapos], resn)
				return {}
			rtdict[k] =  (midict[mikey], resn)

		return rtdict


	# get corresponding map from sequence postion to msa position
	# idxMap[seq_index] = msa_index
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

		# msai2seqi: 3128:244
		msai2seqi = {}
		for k in idxMap:
			if idxMap[k] in msai2seqi:
				print 'getPosMap()::Error!: duplicate msa index: %d' % idxMap[k]
			msai2seqi[idxMap[k]] = k

		seqi2msai = idxMap

		return (seqi2msai, msai2seqi)

	def setTarget(self, t):
		for s in self.msaArray:
			if t in s[0]:
				print 'target found: %s' % s[0]
				self.target = s
				return
		print 'target %s not found!' % t
		sys.exit(-1)

	# generate concised msa scoreboard
	# cutoff: 0% ~ 100%, critera for dropping msa columns. how many gaps in the column
	# weight_cutoff: 0 ~ 1.0 hamming distance between two sequences in msa. How much similar between two sequences in msa.
	# 				 within this distance. two sequence are considered to be in a group with the same weight for the frequency calculation
	def msaboard(self, cutoff): #, weight_cutoff):
		print 'Converting msa to scoreboard ...'
		scoreboard = []
		addboard = np.zeros(self.seqlen)
		for s in self.msaArray:
			scoreboard.append([self.scoreValue[a.upper()] for a in s[1]])
			addboard+=np.array([self.scoreBinary[a.upper()] for a in s[1]]) # calculate gap proportion

		# too slow to do it in python
		#self.calcWeight(scoreboard, weight_cutoff)

		# get conserved columns
		if self.target[0] == '':
			print 'target is empty. Cannot reduce columns'
			sys.exit(-1)

		indices = [i for i in xrange(0, self.seqlen) if (addboard[i]/self.seqNum  > cutoff and self.scoreBinary[self.target[1][i].upper()] != 0)]
		#print np.array(scoreboard)
		return (np.array(scoreboard)[:,indices], indices)


	'''
	W = (1./(1+sum(squareform(pdist(encoded_focus_alignment, 'hamm')<theta))));
	'''
	def calcWeight(self, scoreboard, theta):
		print 'Calculating sequence weight ...'

		pdist = np.ones((self.seqNum, self.seqNum)) # count itself out
		for i in xrange(0, self.seqNum):
			for j in xrange(i+1, self.seqNum):
				pdist[i,j] = (np.array(scoreboard[i])!= np.array(scoreboard[j])).mean() # hamming distance
				pdist[j,i] = pdist[i,j]
				print '%d/%d' % (33629*i+j, 33629*33629)

		#print sum(pdist<theta)
		self.weight = 1.0/(1.0+sum(pdist<theta))
		np.savetxt('1k2p_Pf07714_weight.txt', self.weight)


	# write scoreboard
	def writeScoreboard(self, outfile, cutoff, weight_cutoff):
		fout = open(outfile, 'w')
		score, varlist = self.msaboard(cutoff, weight_cutoff)
		fout.write('%s\n' % repr(varlist))
		for line in score:
			fout.write('%s\n' % ','.join([str(i) for i in line]))
		fout.close()


	# reduce a full alignment to a non-redundant alignment
	# all the sequences in the MSA has a hamming distance >= cutoff
	def hammingReduction(self, outfile, cutoff):
		# get non gap positions
		non_gap_pos = []
		target_seq = self.target[1]
		for i in xrange(0, len(target_seq)):
			if self.scoreValue[target_seq[i]] != 0:
				non_gap_pos.append(i)

		print repr(non_gap_pos)
		non_gap_seq = np.array(list(target_seq))[non_gap_pos]
		print ''.join(non_gap_seq)

		nrArray = []
		nrArray.append(self.target)
		# selecting sequences from msaArray to nrArray
		for s in self.msaArray:
			s_seq = np.array(list(s[1]))
			add_flag = True
			for t in nrArray:
				t_seq = np.array(list(t[1]))
				if (s_seq[non_gap_pos]!=t_seq[non_gap_pos]).mean() < 0.38:
					add_flag = False
					break
			if add_flag == True:
				nrArray.append(s)

		print '\nreduced msa %d/%d' % (len(nrArray), len(self.msaArray))
		print 'writing output to [%s] ...' % outfile

		fout = open(outfile, 'w')
		for fa in nrArray:
			fout.write('>%s\n%s\n' % (fa[0], fa[1]))
		fout.close()
		print 'done.'