import sys
import math
import commp as cp
import numpy as np
from msa import msa
from protein import protein

def contact(x1,y1,z1,x2,y2,z2,cutoff=6.5):
	dist =  np.linalg.norm(np.array((x1,y1,z1))-np.array((x2,y2,z2)))
	if dist <= cutoff:
		return True
	else:
		return False

# only work for coarse-grained pdb
# CA, SGC ..
def getResDict(p):
	return dict((a.resSeq, a) for a in p.atoms)

# key: column index
# value: MSA column in AA alphabet
def getColDict(mscore, pmap):
	return dict((v, [cp.scorerdict[a] for a in mscore[:,v]]) for v in pmap.keys())

# v1 v2 : position pair from miplist
# mscore: msascore
# rdict: protein resi -> atom diction
# pmap: map between pdb resi and msa position
# mip : SDII
def report_hatvar(v1,v2, cdict, rdict, pmap, ensemble):
	# check contact
	a1 = rdict[pmap[v1]]
	a2 = rdict[pmap[v2]]
	isContact = contact(a1.x,a1.y,a1.z,a2.x,a2.y,a2.z)

	# column list in AA alphabet
	#c1 =[cp.scorerdict[a] for a in mscore[:,v1]]
	c1 = cdict[v1]
	ht1 = np.array([cp.aaprop[c][21] for c in c1])
	std1 = np.std(ht1)
	#print 'c1: %d\n' % len(c1)

	#c2 =[cp.scorerdict[a] for a in mscore[:,v2]]
	c2 = cdict[v2]
	ht2 = np.array([cp.aaprop[c][21] for c in c2])
	std2 = np.std(ht2)
	#print 'c2: %d\n' % len(c2)

	ht_12 = np.array([(cp.aaprop[c1[i]][21] + cp.aaprop[c2[i]][21]) for i in xrange(0, len(c1))])
	std12 = np.std(ht_12)

	mip = ensemble[0]
	mi = ensemble[1]
	# 13 column
	return '%d %d %d %d %r %.6f %.6f %.4f %.4f %.4f %.4f %.4f %.4f\n' % (v1, pmap[v1], v2, pmap[v2], isContact, mip, mi, std1, std2, std12, (std1+std2-std12), (std1+std2-std12)*(1-mip), (std1+std2-std12)*(1-mi) )


# v1 v2 : position pair from miplist
# mscore: msascore
# rdict: protein resi -> atom diction
# pmap: map between pdb resi and msa position
# mip : SDII
def report_hat_entropy(v1,v2, cdict, rdict, pmap, ensemble):
	# check contact
	a1 = rdict[pmap[v1]]
	a2 = rdict[pmap[v2]]
	isContact = contact(a1.x,a1.y,a1.z,a2.x,a2.y,a2.z)

	# column list in AA alphabet
	#c1 =[cp.scorerdict[a] for a in mscore[:,v1]]
	c1 = cdict[v1]
	ht1 = np.array([cp.aaprop[c][21] for c in c1])
	h1 = cp.entropy([ht1])
	#print 'c1: %d\n' % len(c1)

	#c2 =[cp.scorerdict[a] for a in mscore[:,v2]]
	c2 = cdict[v2]
	ht2 = np.array([cp.aaprop[c][21] for c in c2])
	h2 = cp.entropy([ht2])

	# entropy of total number of heavy atom for both column
	# one column data
	ht_12 = np.array([math.fabs(cp.aaprop[c1[i]][21] - cp.aaprop[c2[i]][21]) for i in xrange(0, len(c1))])
	h12 = cp.entropy([ht_12])

	# the configuration of heavy atoms
	# for example: [4,6] [3,4]
	#hts12 = np.array([float('%s.%s' % ( (cp.aaprop[c1[i]][21], cp.aaprop[c2[i]][21]) if cp.aaprop[c1[i]][21] >= cp.aaprop[c2[i]][21] else (cp.aaprop[c2[i]][21], cp.aaprop[c1[i]][21]) )) for i in xrange(0, len(c1))])
	#print repr(hts12)
	#hs12 = cp.entropy([hts12])
	mip = ensemble[0]
	mi = ensemble[1]
	# 10 column
	return '%d %d %d %d %r %.6f %.6f %.4f %.4f %.4f %.4f\n' % (v1, pmap[v1], v2, pmap[v2], isContact, mip, mi, h1, h2, h12, (h1+h2-h12))



def main():
	if len(sys.argv)<6:
		print '$ python proc_coevol_heavyatom.py P00966.pdb-PF00764.map PF00764_p90.txt P00966_A.sgc.pdb PF00764_p90.mip outfile'
		return 

	mapfile = sys.argv[1]
	msafile = sys.argv[2]
	pdbfile = sys.argv[3]
	mipfile = sys.argv[4]
	outfile = sys.argv[5]

	# load mip
	# 43-44 0.43466620 0.12479602 0.30987018
	print 'loading mip...'
	miplist = []
	with open(mipfile) as fp:
		for line in fp:
			sarr= line.split()
			idxarr = sarr[0].split('-')
			miplist.append((int(idxarr[0]), int(idxarr[1]), float(sarr[3]), float(sarr[1])  ))
	print 'sorting mip ...'

	#miplist.sort(key=lambda tup: tup[2], reverse=True)
	#maxmip = miplist[0][2]

	maxmip = max(miplist,key=lambda x:x[2])[2]
	maxmi = max(miplist,key=lambda x:x[3])[3]
	print 'max mip: %.4f' % maxmip
	print 'max mi: %.4f' % maxmi
	#for t in miplist:
	#	print '%s %s %s' % (t[0], t[1], t[2])

	# load pdb - > resi map
	# 10 A 44 A
	print 'loading msa resi map ...'
	pairmap = []
	with open(mapfile) as fp:
		for line in fp:
			sarr = line.split()
			pairmap.append((int(sarr[0]), int(sarr[2])))
	#for k,v in pairmap:
	#	print k,v
	msa2resi = dict((v,k) for k,v in pairmap)

	print 'load msa ...'
	m = msa(msafile)
	score = m.msascore()
	cdict = getColDict(score, msa2resi)
	print score.shape

	print 'loading protein ...'
	sgc = protein(pdbfile)
	rdict = getResDict(sgc)
	#print repr(getResDict(sgc))

	print 'generating report ...'
	output = ''.join([report_hatvar(mip[0],mip[1], cdict, rdict, msa2resi, [mip[2]/maxmip,mip[3]/maxmi]) for mip in miplist if (mip[0] in msa2resi) and (mip[1] in msa2resi)])
	#output = ''.join([report_hat_entropy(mip[0],mip[1], cdict, rdict, msa2resi, [mip[2]/maxmip,mip[3]/maxmi]) for mip in miplist if (mip[0] in msa2resi) and (mip[1] in msa2resi)])
	with open(outfile, 'w') as fp:
		fp.write(output)
	print 'done. save to %s' % outfile

if __name__ == '__main__':
	main()