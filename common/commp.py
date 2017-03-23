import sys
import math

# given two lists of coordinates. {1,..,i,..., n} in [x,y,z] format
# return RMSD
# RMSD = sqrt( 1/n * \sum_i (|| v_i - w_i ||^2)   )
def rmsd(v, w):
	if len(v) != len(w):
		print 'error: vector length mismatch. v: %d w: %d' % (len(v), len(w))
		exit(1)
	print repr(v), repr(w)
	d = [((v[i][0]-w[i][0])*(v[i][0]-w[i][0]) + (v[i][1]-w[i][1])*(v[i][1]-w[i][1]) + (v[i][2]-w[i][2])*(v[i][2]-w[i][2])) for i in xrange(0, len(v))] 
	#print repr(d)
	return math.sqrt(sum(d)/len(d))


# given two strings
# normal sequence & aligned sequence
# return map 1. key=1  pos[s1] = s2; 2. key=2 pos[s2] = s1
# s1: aligned string index, s2: pdb sequence index
def posmap(s1, s2, key=1):
	gap = ['.', '-', '_']
	ps1 = s1.translate(None, ''.join(gap))
	ps2 = s2.translate(None, ''.join(gap))
	#print 'ps1: %s\nps2: %s' % (ps1, ps2)

	if ps1!=ps2:
		print 'error: not homo-str'
		print 'ps1: %s\nps2: %s' % (ps1, ps2)
		exit(1)

	i=0
	j=0
	retmap={}
	while(i<len(s1) and j<len(s2)):
		if s1[i] in gap:
			i+=1
			continue
		if s2[j] in gap:
			j+=1
			continue
		if s1[i]==s2[j]:
			if key == 1:
				retmap[i] = j
			else:
				retmap[j] = i
			i+=1
			j+=1

	if len(retmap)!=len(ps1):
		print 'error: incomplete map: len:%d, ps len: %d' % (len(retmap), len(ps1))
		exit(1)

	'''
	print 's1:%s\ns2:%s' % (s1,s2)
	for k in retmap:
		if key == 1:
			print 's1[%d]:%s = s2[%d]:%s' % (k, s1[k], retmap[k], s2[retmap[k]])
		else:
			print 's2[%d]:%s = s1[%d]:%s' % (k, s2[k], retmap[k], s1[retmap[k]])
	'''
	return retmap
