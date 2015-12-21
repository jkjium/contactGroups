# contact group class
# 2ztc_C.domain,LAAAL,155 188 165 168 192,5

class cgroup(object):
	def __init__(self, strline=''):
		if strline == '':
			self.pdb = ''
			self.AAgroup = ''
			self.resi=[]
		else:
			strArray = strline.split(',')
			self.pdb = strArray[0]
			self.AAgroup = strArray[1]
			self.resi = strArray[2].split(' ')

		self.alphabet = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
		self.scoreboard = {'A':0,'C':0,'D':0,'E':0,'F':0,'G':0,'H':0,'I':0,'K':0,'L':0,'M':0,'N':0,'P':0,'Q':0,'R':0,'S':0,'T':0,'V':0,'W':0,'Y':0}
		self.fillScoreboard()


	def getString(self):
		return ('%s,' + 
		 	  '%s,' +
		 	  '%s\n') % \
		 	  (
		 	  	self.pdb,
		 	  	self.AAgroup,
		 	  	' '.join(self.resi)
		 	  )   

	def getSize(self):
		return len(self.AAgroup)

	def fillScoreboard(self):
		for A in self.AAgroup:
			self.scoreboard[A]+=1

	def scoreboard2str(self):
		Alist = []
		for A in self.alphabet:
			Alist.append(str(self.scoreboard[A]))
		return ','.join(Alist)

