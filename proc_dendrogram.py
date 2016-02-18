import sys
import numpy as np
from scipy.cluster.hierarchy import dendrogram
from scipy.cluster.hierarchy import linkage

from atom import atom
from AAmap import AAmap
from protein import protein
#import matplotlib.pyplot as plt

# contact group parser
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

	def writeLeaves(self, resimap):
		resi = []
		resn = []
		for i in self.leaves:
			resi.append(resimap[i][0])
			resn.append(resimap[i][1])
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


# visualization function
# draw cluster linkage
def augmented_dendrogram(*args, **kwargs):
	ddata=dendrogram(*args, **kwargs)
	if not kwargs.get('no_plot', False):
		for i,d in zip(ddata['icoord'], ddata['dcoord']):
			x=0.5*sum(i[1:3])
			y=d[1]
			plt.plot(x,y,'ro')
			str = '%.3g' % y
			plt.annotate(str, (x,y), xytext=(0,-8), textcoords='offset points', va='top', ha='center')
	return ddata


def checkProximity(h, pdist, cutoff):
	#print 'check cluster %d' % (h.clusterID)
	for i in xrange(0, len(h.leaves)):
		for j in xrange(i+1, len(h.leaves)):
			if h.leaves[i] < h.leaves[j]:
				key = '%d-%d' % (h.leaves[i], h.leaves[j])
			elif h.leaves[i] > h.leaves[j]:
				key = '%d-%d' % (h.leaves[j], h.leaves[i])
			elif h.leaves[i] == h.leaves[j]:
				continue
			if pdist[key] > cutoff:
				return False
	return True


def checkProximity2(h1, h2, pdist, cutoff):
	#print 'check c1: %d, c2: %d' % (h1.clusterID, h2.clusterID)
	for i in h1.leaves:
		for j in h2.leaves:
			if i < j:
				key = '%d-%d' % (i,j)
			elif i > j:
				key = '%d-%d' % (j,i)
			elif i == j:
				continue
			if pdist[key] > cutoff:
				return False
	return True


def main():

	if len(sys.argv) < 3:
		print 'python proc_dendrogram.py preffix cutoff'
		exit 

	preffix = sys.argv[1]
	cutoff = float(sys.argv[2])
	# load tip pdb file
	pr = protein(preffix)
	aamap = AAmap()
	n = len(pr.atoms)

	resimap = {}
	print 'writing %s.resimap ...' % (preffix)
	fr = open(preffix+'.resimap', 'w')
	px = []

	count = 0
	for a in pr.atoms:
		px.append((a.x, a.y, a.z))
		resimap[count] = ('%s%d' % (a.chainID, a.resSeq), aamap.getAAmap(a.resName))
		fr.write('%d %s%d %s\n' % (count, a.chainID, a.resSeq, aamap.getAAmap(a.resName)))
		count+=1
	fr.close()

	x = np.array(px)

	# calculate pairwised distance
	pdist = {}
	print 'writing %s.pdist ...' % (preffix)
	fo=open(preffix+'.pdist','w')
	for i in xrange(0,len(x)):
		for j in xrange(i+1,len(x)):
			dist = np.linalg.norm(x[i]-x[j])
			pdist['%d-%d' % (i,j)] = dist
			fo.write('%d-%d : %f\n' % (i,j,dist))
	fo.close()

	# for hc extraction
	hcdict = {}
	hclist = []
	existdict = {}

	#linkage_matrix = linkage(x, "single")
	linkage_matrix = linkage(x, "complete")
	#ddata = augmented_dendrogram(linkage_matrix, color_threshold=1)
	#plt.show()
	print 'writing %s.hcluster ...' % (preffix)
	fo1 = open(preffix+'.hcluster', 'w')
	m = linkage_matrix
	for i in xrange(0,len(m)):
		#print '%d %d %d %f %d' % (n+i,m[i,0],m[i,1],m[i,2],m[i,3])
		fo1.write('%d %d %d %f %d\n' % (n+i,m[i,0],m[i,1],m[i,2],m[i,3]))
		hcline = '%d %d %d %f %d' % (n+i,m[i,0],m[i,1],m[i,2],m[i,3])
		h = hc(hcline, n)
		hcdict[h.clusterID] = h
		hclist.append(h)		
	fo1.close()

	# resolve leaves for each cluster
	print 'resolving leaves ...'
	for h in hclist:
		h.getChildren(hcdict)
		#h.dump()


	print 'iterating clusters for largest proximity contact ...'
	for i in xrange(0, n):
		leafstr = '%d %d %d 0.0 1' % (i, i, i)
		h = hc(leafstr, n)
		h.leaves = [i]
		hcdict[i] = h
		#hcdict[i].dump()


	# add single leaf in
	for i in xrange(0, n):
		existdict[i]= True

	for h in hclist:
		if h.dist <= cutoff:
			if h.c1 in existdict and h.c2 in existdict: # both been checked before
				#print '1AA'
				if existdict[h.c1] == True and existdict[h.c2] == True:
					ret = checkProximity2(hcdict[h.c1], hcdict[h.c2], pdist, cutoff)
					existdict[h.clusterID] = ret
					if ret == True: # combine both and delete sub cluster in the dict
						existdict[h.c1] = False
						existdict[h.c2] = False
				elif existdict[h.c1] == False or existdict[h.c2] == False:
					existdict[h.clusterID] = False

			elif h.c1 in existdict and h.c2 not in existdict:
				#print '1AB'
				if existdict[h.c1] == False: # c1 is not a contact; get h
					existdict[h.clusterID] = False
					existdict[h.c2] = checkProximity(hcdict[h.c2], pdist, cutoff) # get c2
				elif existdict[h.c1] == True: # c1 is a contact; get c2 then get h = c1 and c2
					ret = checkProximity(hcdict[h.c2], pdist, cutoff) # get c2
					existdict[h.c2] = ret
					if ret == False:
						existdict[h.clusterID] = False
					elif ret == True: # h.c2 is a contact
						ret1 = checkProximity2(hcdict[h.c1], hcdict[h.c2], pdist, cutoff)
						existdict[h.clusterID] = ret1
						if ret1 == True:
							existdict[h.c1] = False
							existdict[h.c2] = False

			elif h.c1 not in existdict and h.c2 in existdict:
				#print '1BA'
				if existdict[h.c2] == False: # c2 is not a contact; get h
					existdict[h.clusterID] = False
					existdict[h.c1] = checkProximity(hcdict[h.c1], pdist, cutoff) # get c1
				elif existdict[h.c2] == True: # c2 is a contact; get c1 then get h = c1 and c2
					ret = checkProximity(hcdict[h.c1], pdist, cutoff) # get c1
					existdict[h.c1] = ret
					if ret == False:
						existdict[h.clusterID] = False
					elif ret == True: # h.c1 is a contact
						ret1 = checkProximity2(hcdict[h.c1], hcdict[h.c2], pdist, cutoff)
						existdict[h.clusterID] = ret1
						if ret1 == True:
							existdict[h.c1] = False
							existdict[h.c2] = False

			elif h.c1 not in existdict and h.c2 not in existdict:
				#print '1BB'
				r1 = checkProximity(hcdict[h.c1], pdist, cutoff)
				existdict[h.c1] = r1
				r2 = checkProximity(hcdict[h.c2], pdist, cutoff)
				existdict[h.c2] = r2
				if r1 == False or r2 == False:
					existdict[h.clusterID] = False
				elif r1 == True and r2 == True:
					ret = checkProximity2(hcdict[h.c1], hcdict[h.c2], pdist, cutoff)
					if ret == True:
						existdict[h.c1] = False
						existdict[h.c2] = False

		elif h.dist > cutoff:
			#print '0XX'
			existdict[h.clusterID] = False
			if h.c1 not in existdict:
				existdict[h.c1] = checkProximity(hcdict[h.c1], pdist, cutoff)
			if h.c2 not in existdict:
				existdict[h.c2] = checkProximity(hcdict[h.c2], pdist, cutoff)


	# print out the result
	print 'writing result into %s.hcg' % preffix  
	fout = open(preffix+'.hcg', 'w')
	count=0
	for hid in existdict:
		#if hid >= N and existdict[hid] == True:
		if existdict[hid] == True:
			#fout.write('%d: %r, %s' % (hid, existdict[hid], hcdict[hid].writeString()))
			fout.write('%s,%s\n' % (preffix, hcdict[hid].writeLeaves(resimap)))
			count+=len(hcdict[hid].leaves)
	print '%d leaves in total\n' % count

if __name__ == '__main__':
	main()
