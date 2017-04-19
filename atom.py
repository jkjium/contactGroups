
__all__=['atom']

class atom(object):


    def __init__(self, line):
        line=line.ljust(81).strip()
        self.inputStr=line
        self.outputStr=""
        self.recName=line[0:6]
		#self.recName="ATOM  "
        self.serial=int(line[6:11])
        self.name=line[12:16]
        self.altLoc=line[16]
        self.resName=line[17:20]
        self.chainID=line[21]
        self.resSeq=int(line[22:26])
        self.iCode=line[26]
        self.x=float(line[30:38])
        self.y=float(line[38:46])
        self.z=float(line[46:54])
        self.occupancy=float(line[54:60])
        self.tempFactor=float(line[60:66])
        self.element=line[76:78]
        self.charge=line[78:]
  
    def getCoor(self):
        return [self.x,self.y,self.z]


    def dump(self):
         outStr=self.writeAtom()
         print '%s' % outStr,
		#print self.inputStr	
		#print self.outputStr	
		# print ('recName:[%s]\n' + 
		# 	  'serial:[%d]\n' +
		# 	  'name:[%s]\n' +
		# 	  'altLoc:[%s]\n' +
		# 	  'resName:[%s]\n' +
		# 	  'chainID:[%s]\n' +
		# 	  'resSeq:[%d]\n' +
		# 	  'iCode:[%s]\n' +
		# 	  'x:[%f]\n' +
		# 	  'y:[%f]\n'+
		# 	  'z:[%f]\n' +
		# 	  'occupancy:[%f]\n' +
		# 	  'tempFactor:[%f]\n' +
		# 	  'element:[%s]\n' +
		# 	  'charge:[%s]\n') % \
		# 	  (
		# 	  	self.recName,
		# 	  	self.serial,
		# 	  	self.name,
		# 	  	self.altLoc,
		# 	  	self.resName,
		# 	  	self.chainID,
		# 	  	self.resSeq,
		# 	  	self.iCode,
		# 	  	self.x,
		# 	  	self.y,
		# 	  	self.z,
		# 	  	self.occupancy,
		# 	  	self.tempFactor,
		# 	  	self.element,
		# 	  	self.charge
		# 	  )

    def writeAtom(self):
		self.outputStr = ('%s%s %s%s%s %s%s%s   %s%s%s%s%s          %s%s\n') % \
			(
			  	self.recName,
			  	str(self.serial).rjust(5),
			  	self.name,
			  	self.altLoc,
			  	self.resName,
			  	self.chainID,
			  	str(self.resSeq).rjust(4),
			  	self.iCode,
			  	str('%.3f' % (self.x)).rjust(8),
			  	str('%.3f' % (self.y)).rjust(8),
			  	str('%.3f' % (self.z)).rjust(8),
			  	str('%.2f' % (self.occupancy)).rjust(6),
			  	str('%.2f' % (self.tempFactor)).rjust(6),
			  	self.element,
			  	self.charge
			  )
		return self.outputStr

