# contact group class
# 2ztc_C.domain,LAAAL,155 188 165 168 192,5
# 1t3r.tip,LTT,B97 A26 B26
import operator

class cg(object):
	def __init__(self, strline, var_alphabet):
		strArray = strline.split(',')
		self.pdb = strArray[0]
		self.AAgroup = strArray[1]
		self.resi = strArray[2].split(' ')
		self.AAdict = {}
		for i in xrange(0,len(self.AAgroup)):
			# self.AAdict['B97'] = 'L'
			self.AAdict[self.resi[i]] = self.AAgroup[i]

		self.sortedAA = []
		self.sortedResi = []
		self.sortAAgroup()
		#print self.sortedAA, self.sortedResi

		# for AA type: H + - P 'C,G,O'
		# O for proline
		self.AAType = []

		self.scoreboard = {}
		self.alphabet = var_alphabet
		for a in self.alphabet:
			self.scoreboard[a] = 0

		#self.alphabet = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
		#self.scoreboard = {'A':0,'C':0,'D':0,'E':0,'F':0,'G':0,'H':0,'I':0,'K':0,'L':0,'M':0,'N':0,'P':0,'Q':0,'R':0,'S':0,'T':0,'V':0,'W':0,'Y':0}

	# output scoreboard for AA + accesibility type contact group
	def nascore(self, na):
		# key = '%s%s%s' % (aamap.getAAmap(r.resn), r.chain, r.resi)
		for k in self.AAdict:
			residue = self.AAdict[k]
 			na_key = '%s%s' % (residue, k)
 			varname = '%s%s' % (residue, na.accessible(na_key))
 			self.scoreboard[varname]+=1

		retlist = []
		for a in self.alphabet:
			retlist.append(str(self.scoreboard[a]))
		return ','.join(retlist)


	# put sorted cg type string
	def getType(self, typeMap):
		for A in self.sortedAA:
			self.AAType.append(typeMap[A])
		self.AAType.sort()


	# put CG and indices in sorted Alphabet
	def sortAAgroup(self):
		sorted_dict = sorted(self.AAdict.items(), key=operator.itemgetter(1))
		for pair in sorted_dict:
			self.sortedAA.append(pair[1])
			self.sortedResi.append(pair[0])



	def getString(self):
		return ('%s,' + 
		 	  #'%s,' +
		 	  #'%s,' +
		 	  '%s,' +
		 	  '%s,' +
		 	  '%s\n') % \
		 	  (
		 	  	self.pdb,
		 	  	#self.AAgroup,
		 	  	#' '.join(self.resi),
		 	  	''.join(self.sortedAA),
		 	  	' '.join(self.sortedResi),
		 	  	''.join(self.AAType)
		 	  )   

	def getSize(self):
		return len(self.AAgroup)

	def fillScoreboard(self):
		for A in self.AAgroup:
			self.scoreboard[A]+=1

	def scoreboard2str(self):
		self.fillScoreboard()
		Alist = []
		for A in self.alphabet:
			Alist.append(str(self.scoreboard[A]))
		return ','.join(Alist)

