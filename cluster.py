from atom import atom
from AAmap import AAmap


__all__=['cluster']

class cluster(object):

    def __init__(self, pdb, top, pfam, str, pdbidx, seqheader, alignstr, alignidx, center, cutoff, scutoff, flag, weight, desc):
        self.pdb=pdb
        self.top=top
        self.pfam=pfam

        self.str=str
        self.pdbidx=pdbidx
   
        self.seqheader=seqheader
        self.alignstr=alignstr
        self.alignidx=alignidx

        self.center=center
        self.cutoff=cutoff
        self.scutoff=scutoff

        self.flag=flag
        self.weight = weight
        self.desc=desc
        
        self.typeStr=''
        self.pdbResSeq=''
        self.thetaPhi = []
        
    def getString(self):
 		return ('%s_' + 
		 	  '%s_' +
		 	  '%s_' +
		 	  '%s_' +
		 	  '%s_' +
		 	  '%s_' +
		 	  '%s_' +
		 	  '%s_' +
		 	  '%s_' +
		 	  '%f_' +
		 	  '%d_' +
		 	  '%d_' +
                 '%f_' +
		 	  '%s\n') % \
		 	  (
		 	  	self.pdb,
		 	  	self.top,
		 	  	self.pfam,
		 	  	self.str,
		 	  	self.pdbidx,
                    self.seqheader,
		 	  	self.alignstr,
		 	  	self.alignidx,
		 	  	self.center,       
		 	  	self.cutoff,
		 	  	self.scutoff,
		 	  	self.flag,
                    self.weight,
		 	  	self.desc
		 	  )   


    def addNeighbor(self, amap, at, index):
        self.str=('%s%s') % (self.str, amap.getAAmap(at.resName))
        self.typeStr = ('%s%s') % (self.typeStr, amap.getA2Tmap(amap.getAAmap(at.resName)))
        self.pdbidx = ('%s %s') % (self.pdbidx, str(index))
        self.pdbResSeq = ('%s %d') % (self.pdbResSeq, at.resSeq)
    
    def getPDBidxArray(self):
        s=self.pdbidx[1:len(self.pdbidx)]
        return s.split(' ')

    def toString(self):
 		return ('%s|' + 
		 	  '%s|' +
		 	  '%s|' +
		 	  '%s|' +
		 	  '%s|' +
		 	  '%s|' +
		 	  '%s|' +
		 	  '%f|' +
		 	  '%d|' +
		 	  '%s|'+
		 	  '%s|' +
		 	  '%d|' +
                 '%d|' +
                 '%f|' +
		 	  '%s|\n') % \
		 	  (
		 	  	self.pdb,
		 	  	self.top,
		 	  	self.pfam,
		 	  	self.str,
		 	  	self.pdbidx,
		 	  	self.alignstr,
		 	  	self.alignidx,
		 	  	self.cutoff,
		 	  	self.scutoff,
		 	  	self.center,
		 	  	self.seqheader,
		 	  	self.flag,
                    len(self.str),
                    self.weight,
		 	  	self.desc
		 	  )      

    def toDBString(self):
 		return ('(\'%s\',' + 
		 	  '\'%s\',' +
		 	  '\'%s\',' +
		 	  '\'%s\',' +
		 	  '\'%s\',' +
		 	  '\'%s\',' +
		 	  '\'%s\',' +
		 	  '\'%s\','+
		 	  #'\'%s\',' +      
		 	  #'%f,' +
		 	  #'%d,' +
		 	  '%d,' +
                 '%d,' +
                 '%f,' +
		 	  '\'%s\')') % \
		 	  (
		 	  	self.pdb,
		 	  	self.top,
		 	  	self.pfam,
		 	  	self.str,
		 	  	self.pdbidx,
		 	  	self.seqheader,       
		 	  	self.alignstr,
		 	  	self.alignidx,
		 	  	#self.center,       
		 	  	#self.cutoff,
		 	  	#self.scutoff,
		 	  	self.flag,
                    len(self.str),
                    self.weight,
		 	  	self.desc
		 	  ) 
       
    def dump(self):
 		print ('pdb:[%s]\n' + 
		 	  'top:[%s]\n' +
		 	  'pfam:[%s]\n' +
		 	  'str:[%s]\n' +
		 	  'pdbidx:[%s]\n' +
		 	  'alignstr:[%s]\n' +
		 	  'alignidx:[%s]\n' +
		 	  'cutoff:[%f]\n' +
		 	  'scutoff:[%d]\n' +
		 	  'center:[%s]\n'+
		 	  'seqheader:[%s]\n' +
		 	  'flag:[%d]\n' +
                 'len:[%d]\n' +
		 	  'desc:[%s]\n') % \
		 	  (
		 	  	self.pdb,
		 	  	self.top,
		 	  	self.pfam,
		 	  	self.str,
		 	  	self.pdbidx,
		 	  	self.alignstr,
		 	  	self.alignidx,
		 	  	self.cutoff,
		 	  	self.scutoff,
		 	  	self.center,
		 	  	self.seqheader,
		 	  	self.flag,
		 	  	self.desc
		 	  )       
    #def init_by_string(self, line):
        