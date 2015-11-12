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

# get x
fp = open('1oai.coor', 'r')
px = []
lines = fp.readlines()
n=len(lines)
for line in lines:
	strArr = line.strip().split(' ')
	px.append((float(strArr[0]), float(strArr[1]), float(strArr[2])))
x = np.array(px)
fp.close()

# calculate pairwised distance
print('writing pairwised distance file ...')
fo=open('pdist.txt','w')
for i in xrange(0,len(x)):
	for j in xrange(i+1,len(x)):
		dist = np.linalg.norm(x[i]-x[j])
		fo.write('%d-%d : %f\n' % (i,j,dist))
fo.close()

linkage_matrix = linkage(x, "single")
#ddata = augmented_dendrogram(linkage_matrix, color_threshold=1)
#plt.show()
print('writing hcluster file ...')
fo1 = open('hcluster.txt', 'w')
m = linkage_matrix
for i in xrange(0,len(m)):
	#print '%d %d %d %f %d' % (n+i,m[i,0],m[i,1],m[i,2],m[i,3])
	fo1.write('%d %d %d %f %d\n' % (n+i,m[i,0],m[i,1],m[i,2],m[i,3]))
fo1.close()
