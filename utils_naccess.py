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
			# VB1039
			key1 = '%s%s%s' % (sarr[4], sarr[1], sarr[2])
			key2 = '%s%s%s' % (sarr[5], sarr[1], sarr[3])

			rel1 = na.rsaDict[key1].AA_REL if key1 in na.rsaDict else '-191'
			rel2 = na.rsaDict[key2].AA_REL if key2 in na.rsaDict else '-191'

			outlist.append('%s %s %s' % (line, rel1, rel2))

	with open(outfile, 'w') as fout:
		fout.write('%s' % '\n'.join(outlist))
	cp._info('save rcflat %s' % outfile)



def foo(arglist):
	if len(arglist) < 1:
		cp._err('Usage: python utils_naccess.py foo rsafile')
	na = naccess(arglist[0])
	na.dump()

if __name__ == '__main__':
	cp.dispatch(__name__)