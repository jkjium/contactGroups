'''
get resid list of varname
'''
import sys
from naccess import naccess
from naccess import rsa
import commp as cp

def rcflat(arglist):
	if len(arglist) < 2:
		cp._err('Usage: python utils_naccess.py rcflat rsafile cflatfile')

	rsafile = arglist[0]
	cflatfile = arglist[1]
	outfile = cflatfile[:-4] + '.rcflat'

	outlist = []
	na = naccess(rsafile)

	with open(cflatfile) as fp:
		for line in fp:
			# p.pdb,chainid,r1,r2,res1,res2,dist_sgc,dist_tip,dist_ca,pfamid,p1,p2 mip dca rel1 rel2
			# 3h0g B 30 31 L A 5.90678710 6.50986125 3.78331297 PF04563 174 175 -191 56.11176093
			line = line.strip()
			if len(line) == 0:
				continue
			sarr = line.split()
			# LB30
			key1 = '%s%s%s' % (sarr[4], sarr[1], sarr[2])
			key2 = '%s%s%s' % (sarr[5], sarr[1], sarr[3])

			rel1 = na.rsaDict[key1].AA_REL if key1 in na.rsaDict else '-191'
			rel2 = na.rsaDict[key2].AA_REL if key2 in na.rsaDict else '-191'

			outlist.append('%s %s %s' % (line, rel1, rel2))

	with open(outfile, 'w') as fout:
		fout.write('%s\n' % '\n'.join(outlist))
	cp._info('save rcflat %s' % outfile)

# input resi columns
# output surface accessible area % from NAccess
# naccess key: = '%s%s%s' % (aamap.getAAmap(r.resn), r.chain, r.resi)
def rv62rsa(args):
	assert len(args) == 3, 'Usage: python utils_naccess.py res2rsa rsafile chain_resi_resn.vec6 outfile'

	rsafile = args[0]
	resifile = args[1]
	outfile = args[2]

	na = naccess(rsafile)
	outlist = []
	# t: 0.chain1 1.chain2 2.resi1 3.resi2 4.resn1 5.resn2
	for t in cp.loadtuples(resifile):
		k1 = '%s%s%s' % (t[4],t[0],t[2])
		k2 = '%s%s%s' % (t[5],t[1],t[3])
		outlist.append('%.2f %.2f' % (na.rsaDict[k1].AA_REL, na.rsaDict[k2].AA_REL))

	with open(outfile, 'w') as fout:
		fout.write('%s\n' % '\n'.join(outlist))
	cp._info('save rsa vec2 to %s' % outfile)


def foo(arglist):
	if len(arglist) < 1:
		cp._err('Usage: python utils_naccess.py foo rsafile')
	na = naccess(arglist[0])
	na.dump()

if __name__ == '__main__':
	cp.dispatch(__name__)