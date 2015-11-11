from hc import hc
hcdict = {}
hclist = []
fp = open('cluster.txt', 'r')
for line in fp.readlines():
	h = hc(line.strip(), 59)
	hcdict[h.clusterID] = h
	hclist.append(h)
fp.close()

for h in hclist:
	h.getChildren(hcdict)
	h.dump()
