import sys
import numpy as np
from scipy.cluster.hierarchy import dendrogram
from scipy.cluster.hierarchy import linkage
from scipy.cluster.hierarchy import fcluster

import matplotlib.pyplot as plt
import commp as cp

def fancy_dendrogram(*args, **kwargs):
    max_d = kwargs.pop('max_d', None)
    if max_d and 'color_threshold' not in kwargs:
        kwargs['color_threshold'] = max_d
    annotate_above = kwargs.pop('annotate_above', 0)

    ddata = dendrogram(*args, **kwargs)

    if not kwargs.get('no_plot', False):
        plt.title('Hierarchical Clustering Dendrogram (truncated)')
        plt.xlabel('sample index or (cluster size)')
        plt.ylabel('distance')
        for i, d, c in zip(ddata['icoord'], ddata['dcoord'], ddata['color_list']):
            x = 0.5 * sum(i[1:3])
            y = d[1]
            if y > annotate_above:
                plt.plot(x, y, 'o', c=c)
                plt.annotate("%.3g" % y, (x, y), xytext=(0, -5),
                             textcoords='offset points',
                             va='top', ha='center')
        if max_d:
            plt.axhline(y=max_d, c='k')
    return ddata

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

# x: str list in np.array format
#$ python proc_hammingweight.py t_hamming_weight.score 0.9
#array([1, 2, 2, 3, 4, 1, 1], dtype=int32)
#defaultdict(<type 'int'>, {1: 3, 2: 2, 3: 1, 4: 1})
# 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 : 0.333
# 1, 1, 1, 1, 1, 1, 1, 1, 3, 3 : 0.500
# 1, 1, 1, 1, 1, 1, 1, 1, 3, 3 : 0.500
# 1, 1, 1, 1, 1, 3, 3, 3, 4, 4 : 1.000
# 1, 1, 1, 1, 1, 3, 3, 3,15,15 : 1.000
# 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 : 0.333
# 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 : 0.333
def hamming_weight(x, max_d):
	linkage_matrix = linkage(x, "single", metric='hamming')
	ddata = augmented_dendrogram(linkage_matrix, color_threshold=1)
	plt.show()	
	clusters = fcluster(linkage_matrix, max_d, criterion='distance')
	#print repr(clusters)
	normdict = cp.freq(clusters)
	#print repr(normdict)
	weight = [(1.0/normdict[k]) for k in clusters]
	return weight

def main():
	if len(sys.argv) < 3:
		cp._err('Usage: python proc_hammingweight.py PF00000.score similarity_value')

	datafile = sys.argv[1]
	svalue = float(sys.argv[2])
	outfile = '%s.%2d.w' % (datafile, svalue*100)

	score = np.loadtxt(datafile, delimiter=',')
	#with open(datafile) as fp:
	#	strlist = np.array([[cp.aascore['aa'][a] for a in line.strip()] for line in fp if len(line.strip())!=0])
	w = hamming_weight(score, 1-svalue)
	print repr(w)
	'''
	for i in xrange(0, len(w)):
		print '%s : %.03f' % (','.join(['%2d' % v for v in score[i,:]]), w[i])
	print outfile
	np.savetxt(outfile, w)
	#with open(outfile ,'w') as fp:
	#	fp.write(','.join(['%.8f' % v for v in w]))
	cp._info('save weight to %s' % outfile)
	'''
if __name__ == '__main__':
	main()