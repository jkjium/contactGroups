import sys
import numpy as np
from scipy.cluster.hierarchy import dendrogram
from scipy.cluster.hierarchy import linkage
#import matplotlib.pyplot as plt


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

if len(sys.argv) < 2:
	print 'python dendrogram.py preffix'
	return 

# get x
fp = open(preffix+'.spdb', 'r')
px = []
lines = fp.readlines()
fp.close()
n=len(lines)
print 'writing %s.resimap ...' % (preffix)
fr = open(preffix+'.resimap', 'w')
count = 0
for line in lines:
	strArr = line.strip().split(' ')
	px.append((float(strArr[0]), float(strArr[1]), float(strArr[2])))
	fr.write('%d %d %s\n' % (count, int(strArr[3]), strArr[4]))
	count+=1
fr.close()

x = np.array(px)

# calculate pairwised distance
print 'writing %s.pdist ...' % (preffix)
fo=open(preffix+'.pdist','w')
for i in xrange(0,len(x)):
	for j in xrange(i+1,len(x)):
		dist = np.linalg.norm(x[i]-x[j])
		fo.write('%d-%d : %f\n' % (i,j,dist))
fo.close()

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
fo1.close()
