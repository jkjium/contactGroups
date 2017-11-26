import sys
import commp as cp
from protein import protein
from atom import atom

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
	cp._info('%s:%d lines loaded.' % (mapfile, len(resimap)))

	p = protein(pdbfile, chain=chainid)
	sgc = dict((a.resSeq, a) for a in p.atomsbyscgmcenter())
	tip = dict((a.resSeq, a) for a in p.atomsbytip())
	ca = dict((a.resSeq, a) for a in p.ca)
	cp._info('%s-%s:%d sgc %d tip %d ca loaded.' % (pdbfile, chainid, len(sgc), len(tip), len(ca)))

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
				cp._info('err:mismatch residue in %s: %s - sgc[%d]:%s, %s - sgc[%d]:%s' % (mapfile, res1, r1, sgc[r1].resName, res2, r2, sgc[r2].resName))
				return
			dist_sgc = cp.dist([sgc[r1].x, sgc[r1].y, sgc[r1].z], [sgc[r2].x, sgc[r2].y, sgc[r2].z])

			# tip distance
			if res1!=cp.aa2a[tip[r1].resName] or res2!=cp.aa2a[tip[r2].resName]:
				cp._info('err:mismatch residue in %s: %s - tip[%d]:%s, %s - tip[%d]:%s' % (mapfile, res1, r1, tip[r1].resName, res2, r2, tip[r2].resName))
				return
			dist_tip = cp.dist([tip[r1].x, tip[r1].y, tip[r1].z], [tip[r2].x, tip[r2].y, tip[r2].z])

			# ca distance
			if res1!=cp.aa2a[ca[r1].resName] or res2!=cp.aa2a[ca[r2].resName]:
				cp._info('err:mismatch residue in %s: %s - ca[%d]:%s, %s - ca[%d]:%s' % (mapfile, res1, r1, ca[r1].resName, res2, r2, ca[r2].resName))
				return
			dist_ca = cp.dist([ca[r1].x, ca[r1].y, ca[r1].z], [ca[r2].x, ca[r2].y, ca[r2].z])
			#        pdb r1 r2 res1 res2 dist.sgc dist.tip dist.ca pfamid p1 p2
			outstr = '%s %s %d %d %s %s %.8f %.8f %.8f %s %d %d\n' % (p.pdb,chainid,r1,r2,res1,res2,dist_sgc,dist_tip,dist_ca,pfamid,p1,p2)
			fp.write(outstr)
	fp.close()
	cp._info('save to %s' % outfile)



def main():
	if len(sys.argv)<2:
		cp._err('Usage: python utils_flat.py cmd [args ...]')

	dispatch = {
		'flaten':flaten
	}

	if sys.argv[1] in dispatch:
		dispatch[sys.argv[1]](sys.argv[2:])
	else:
		cp._err('invalid cmd: %s' % sys.argv[1])

if __name__ == '__main__':
	main()