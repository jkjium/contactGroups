import sys
from hc import hc

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
	hcdict = {}
	hclist = []
	existdict = {}
	pdist = {}
	resmap = {}

	if len(sys.argv) < 3:
		print 'python proc_hierachialCG.py preffix cutoff'
		print 'need .resimap .pdist .hcluster pre-calculated by dendrogram.py'
		return

	preffix = sys.argv[1]
	cutoff = float(sys.argv[2])
	print 'preffix: [%s]' % preffix
	print 'cutoff: [%f]' % cutoff

	# load index -> residue 
	print 'loading res map ...'
	fr = open(preffix+'.resimap', 'r')
	for line in fr.readlines():
		strArr = line.strip().split(' ')
		resmap[int(strArr[0])] = (int(strArr[1]), strArr[2])
	fr.close()
	N = len(resmap)


	# load pairwised distance
	print 'loading pairwised distance ...'
	fd = open(preffix+'.pdist', 'r')
	for line in fd.readlines():
		strArr = line.strip().split(' ')
		pdist[strArr[0]] = float(strArr[2])
	fd.close()

	# load cluster tree
	print 'loading cluster tree ...'
	fp = open(preffix+'.hcluster', 'r')
	for line in fp.readlines():
		h = hc(line.strip(), N)
		hcdict[h.clusterID] = h
		hclist.append(h)
	fp.close()

	# resolve leaves for each cluster
	print 'resolving leaves ...'
	for h in hclist:
		h.getChildren(hcdict)
		#h.dump()


	print 'iterating clusters for largest proximity contact ...'
	for i in xrange(0, N):
		leafstr = '%d %d %d 0.0 1' % (i, i, i)
		h = hc(leafstr, N)
		h.leaves = [i]
		hcdict[i] = h
		#hcdict[i].dump()


	# add single leaf in
	for i in xrange(0, N):
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
	print 'writing result into %s.cg' % preffix  
	fout = open(preffix+'.hcg', 'w')
	count=0
	for hid in existdict:
		#if hid >= N and existdict[hid] == True:
		if existdict[hid] == True:
			#fout.write('%d: %r, %s' % (hid, existdict[hid], hcdict[hid].writeString()))
			fout.write('%s,%s\n' % (preffix, hcdict[hid].writeLeaves(resmap)))
			count+=len(hcdict[hid].leaves)
	print '%d leaves in total\n' % count



if __name__ == '__main__':
	main()
