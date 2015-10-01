

__all__=['cgroup']

class cgroup(object):
	def __init__(self, strline=''):
		if strline == '':
			self.pdb = ''
			self.AAgroup = ''
			self.resi=[]
		else:
			strArray = strline.split(',')
			self.pdb = strArray[0]
			self.AAgroup = strArray[2]
			self.resi = strArray[4].split(' ')

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

