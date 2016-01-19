# contact group class
# 2ztc_C.domain,LAAAL,155 188 165 168 192,5
import operator

class cgroup(object):
	def __init__(self, strline=''):
		strArray = strline.split(',')
		self.pdb = strArray[0]
		self.AAgroup = strArray[1]
		self.resi = strArray[2].split(' ')
		self.AAdict = {}
		for i in xrange(0,len(self.AAgroup)):
			self.AAdict[self.resi[i]] = self.AAgroup[i]

		self.sortedAA = []
		self.sortedResi = []
		self.sortAAgroup()
		print self.sortedAA, self.sortedResi

		# for AA type: H + - P 'C,G,O'
		# O for proline
		self.AAType = []

		self.alphabet = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
		self.scoreboard = {'A':0,'C':0,'D':0,'E':0,'F':0,'G':0,'H':0,'I':0,'K':0,'L':0,'M':0,'N':0,'P':0,'Q':0,'R':0,'S':0,'T':0,'V':0,'W':0,'Y':0}

	def getType(self, typeMap):
		for A in self.sortedAA:
			self.AAType.append(typeMap[A])


	# put CG and indices in sorted Alphabet
	def sortAAgroup(self):
		sorted_dict = sorted(self.AAdict.items(), key=operator.itemgetter(1))
		for pair in sorted_dict:
			self.sortedAA.append(pair[1])
			self.sortedResi.append(pair[0])



	def getString(self):
		return ('%s\n' + 
		 	  '%s\n' +
		 	  '%s\n' +
		 	  '%s\n' +
		 	  '%s\n' +
		 	  '%s\n') % \
		 	  (
		 	  	self.pdb,
		 	  	self.AAgroup,
		 	  	' '.join(self.resi),
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

