
__all__=['hc']

class hc(object):


    def __init__(self, line, n):

    	self.N = n
    	strArr = line.split(' ')
    	self.clusterID = int(strArr[0])
    	self.c1 = int(strArr[1])
    	self.c2 = int(strArr[2])
    	self.dist = float(strArr[3])
    	self.clusterLen = int(strArr[4])

    	self.leaves = []



    def dump(self):
        outStr=self.writeString()
        print '%s' % outStr

    def writeLeaves(self, resmap):
        resi = []
        resn = []
        for i in self.leaves:
            resi.append(str(resmap[i][0]))
            resn.append(resmap[i][1])
        return '%s,%s' % (''.join(resn), ' '.join(resi))


    def writeString(self):
		return ('%d %d %d %f %d [%d]: %s\n') % (self.clusterID, self.c1, self.c2, self.dist, self.clusterLen, len(self.leaves), str(self.leaves))

    def getLeaves(self, hcdict):
    	ret = []
    	if self.c1 < self.N:
    		ret.append(self.c1)
    	else:
    		ret = ret + hcdict[self.c1].getLeaves(hcdict)

    	if self.c2 < self.N:
    		ret.append(self.c2)
    	else:
    		ret = ret + hcdict[self.c2].getLeaves(hcdict)

    	return ret

    def getChildren(self, hcdict):
    	self.leaves = self.getLeaves(hcdict)
