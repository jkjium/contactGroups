import sys
import commp as cp
from protein import protein
from atom import atom

class cetuple(object):
	def __init__(self, flatstr):
		sarr = flatstr.split(' ')
		#p.pdb,chainid,r1,r2,res1,res2,dist_sgc,dist_tip,dist_ca,pfamid,p1,p2
		self.pdb = sarr[0]
		self.chain = sarr[1]
		self.r1 = sarr[2]
		self.r2 = sarr[3]
		self.rn1 = sarr[4]
		self.rn2 = sarr[5]
		self.cg = {'sgc': float(sarr[6]), 'tip': float(sarr[7]), 'ca': float(sarr[8])}
		self.pfamid = sarr[9]
		self.p1 = sarr[10]
		self.p2 = sarr[11]
		self.ce = [float(e) for e in sarr[12:]] 

	def dump(self):
		cp._info('%s %s %s %s %s %s %s %s %s %s %s' %
			(self.pdb,
			self.chain,
			self.r1,
			self.r2,
			self.rn1,
			self.rn2,
			repr(self.cg),
			self.pfamid,
			self.p1,
			self.p2,
			repr(self.ce))
			)

# python utils_flat.py cecolumn 101m.pdb-A-PF00042.map PF00042_p90.mip 4 std PF00042-std-mip.ce
# cecol starts from 0
def cecolumn(arglist):
	if len(arglist)<5:
		cp._err('Usage: python utils_flat.py cecolumn mapfile cefile cecol method outfile')

	mapfile = arglist[0]
	cefile = arglist[1]
	cecol = int(arglist[2])
	method = arglist[3]
	outfile = arglist[4]

	# load ce values
	cedict = {}
	with open(cefile) as fp:
		for line in fp:
			line = line.strip()
			if len(line) == 0:
				continue
			sarr = line.split(' ')
			k = ('%s %s' % (sarr[0], sarr[1])) if sarr[0] < sarr[1] else ('%s %s' % (sarr[1], sarr[0]))
			v = float(sarr[cecol])
			cedict[k] = v

	# use method to normalize ce values
	if method == 'std':
		normdict = cp.rankstd(cedict)
	else:
		normdict = cedict

	# extract ce for mapped tuples
	resimap = []
	with open(mapfile) as fp:
		for line in fp:
			line = line.strip()
			if len(line) == 0:
				continue
			sarr = line.split(' ')
			#113 L 170 I
			resimap.append(int(sarr[2]))

	fp = open(outfile,'w')
	for n in xrange(0, len(resimap)):
		for m in xrange(n+1, len(resimap)):
			k = '%s %s' % (resimap[n] ,resimap[m])
			outstr = ('%.8f\n' % normdict[k]) if k in normdict else '-191\n'
			#print k,outstr,
			fp.write(outstr)
	fp.close()
	cp._info('save ce to %s' % outfile)


# generate flat stub file
# ec score will be appended
def flaten(arglist):
	if len(arglist) < 5:
		cp._err('Usage: python utils_flat.py mapfile pdbfile chainid PfamId outfile')

	mapfile = arglist[0]
	pdbfile = arglist[1]
	chainid = arglist[2]
	pfamid = arglist[3]
	outfile = arglist[4]

	# load map
	resimap = []
	with open(mapfile) as fp:
		for line in fp:
			line = line.strip()
			if len(line) == 0:
				continue
			sarr = line.split(' ')
			#113 L 170 I
			resimap.append((int(sarr[0]), int(sarr[2]), sarr[1]))
	#cp._info('%s:%d lines loaded.' % (mapfile, len(resimap)))

	p = protein(pdbfile, chain=chainid)
	sgc = dict((a.resSeq, a) for a in p.atomsbyscgmcenter())
	tip = dict((a.resSeq, a) for a in p.atomsbytip())
	ca = dict((a.resSeq, a) for a in p.ca)
	#cp._info('%s-%s:%d sgc %d tip %d ca loaded.' % (pdbfile, chainid, len(sgc), len(tip), len(ca)))

	fp = open(outfile, 'w')
	# pairwise 
	for n in xrange(0, len(resimap)):
		for m in xrange(n+1, len(resimap)):
			i = resimap[n]
			j = resimap[m]

			# first resi
			r1 = i[0]
			p1 = i[1]
			res1 = i[2]

			# second resi
			r2 = j[0]
			p2 = j[1]
			res2 = j[2]

			# sgc distance
			if res1!=cp.aa2a[sgc[r1].resName] or res2!=cp.aa2a[sgc[r2].resName]:
				cp._err('err:mismatch residue in %s: %s - sgc[%d]:%s, %s - sgc[%d]:%s' % (mapfile, res1, r1, sgc[r1].resName, res2, r2, sgc[r2].resName))
			dist_sgc = cp.dist([sgc[r1].x, sgc[r1].y, sgc[r1].z], [sgc[r2].x, sgc[r2].y, sgc[r2].z])

			# tip distance
			if res1!=cp.aa2a[tip[r1].resName] or res2!=cp.aa2a[tip[r2].resName]:
				cp._err('err:mismatch residue in %s: %s - tip[%d]:%s, %s - tip[%d]:%s' % (mapfile, res1, r1, tip[r1].resName, res2, r2, tip[r2].resName))
			dist_tip = cp.dist([tip[r1].x, tip[r1].y, tip[r1].z], [tip[r2].x, tip[r2].y, tip[r2].z])

			# ca distance
			if res1!=cp.aa2a[ca[r1].resName] or res2!=cp.aa2a[ca[r2].resName]:
				cp._err('err:mismatch residue in %s: %s - ca[%d]:%s, %s - ca[%d]:%s' % (mapfile, res1, r1, ca[r1].resName, res2, r2, ca[r2].resName))
			dist_ca = cp.dist([ca[r1].x, ca[r1].y, ca[r1].z], [ca[r2].x, ca[r2].y, ca[r2].z])
			#        pdb r1 r2 res1 res2 dist.sgc dist.tip dist.ca pfamid p1 p2
			outstr = '%s %s %d %d %s %s %.8f %.8f %.8f %s %d %d\n' % (p.pdb,chainid,r1,r2,res1,res2,dist_sgc,dist_tip,dist_ca,pfamid,p1,p2)
			fp.write(outstr)
	fp.close()
	cp._info('save to %s' % outfile)


# extract scol from single condition (cg & ce)
def scolsingle(arglist):
	if len(arglist)< 6:
		cp._err('Usage: python utils_flat.py scolsingle flatfile cgname cgcutoff ceidx cecutoff outfile')

	flatfile = arglist[0]
	cgname = arglist[1] # sgc, tip, ca
	cgcutoff = float(arglist[2])
	ceidx = int(arglist[3])
	cecutoff = float(arglist[4])
	outfile = arglist[5]

	scol = []
	with open(flatfile) as fp:
		for line in fp:
			line = line.strip()
			if len(line)==0:
				continue
			ct = cetuple(line)
			if ct.cg[cgname] <= cgcutoff and ct.ce[ceidx] >= cecutoff:
				#ct.dump()
				scol.append('%s-%s' % (ct.p1, ct.p2))

	with open(outfile, 'w') as fp:
		fp.write(' '.join(scol))
	cp._info('save %d tuples to %s' % (len(scol), outfile))





def main():
	if len(sys.argv)<2:
		cp._err('Usage: python utils_flat.py cmd [args ...]')

	dispatch = {
		'flaten':flaten,
		'cecolumn':cecolumn,
		'scolsingle':scolsingle
	}

	if sys.argv[1] in dispatch:
		dispatch[sys.argv[1]](sys.argv[2:])
	else:
		cp._err('invalid cmd: %s' % sys.argv[1])

if __name__ == '__main__':
	main()