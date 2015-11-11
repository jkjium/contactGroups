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
#np.random.seed(12312)
#n=10
#x=np.random.multivariate_normal([0,0],np.array([[4.0,2.5],[2.5,1.4]]),size=(n,))
# get x end

fp = open('1oai.coor', 'r')
px = []
lines = fp.readlines()
n=len(lines)
for line in lines:
	strArr = line.strip().split(' ')
	px.append((float(strArr[0]), float(strArr[1]), float(strArr[2])))
x = np.array(px)
fp.close()
#print x
#fo=open('dist.txt','w')
#for i in xrange(0,len(x)):
#	for j in xrange(i+1,len(x)):
#		dist = np.linalg.norm(x[i]-x[j])
#		fo.write('%d - %d : %f\n' % (i,j,dist))
#fo.close()
linkage_matrix = linkage(x, "single")
#ddata = augmented_dendrogram(linkage_matrix, color_threshold=1)
#plt.show()
m = linkage_matrix
#print '%f, %f, %f, %f' % (m[0,0],m[0,1],m[0,2],m[0,3])
#print '%f, %f, %f, %f' % (m[1,0],m[1,1],m[1,2],m[1,3])
#print len(m)
for i in xrange(0,len(m)):
	#print (n+i,m[i,0],m[i,1],m[i,2],m[i,3])
	print '%d %d %d %f %d' % (n+i,m[i,0],m[i,1],m[i,2],m[i,3])
