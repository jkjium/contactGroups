'''
20190329 for finding charged buired amino acid pairs (coevolved)
[kjia@lhb-ps1 ~/workspace/pfam31.0/p90_ext] *.rcflat
'''

import commp as cp

def buriedchargepair(arglist):
	if len(arglist) < 2:
		cp._err('Usage: python utils_cflat.py buriedchargepair inrcflatfile outfile')

	infile = arglist[0]
	outfile = arglist[1]

	outlist = []

	# p.pdb,chainid,r1,r2,res1,res2,dist_sgc,dist_tip,dist_ca,pfamid,p1,p2,mip,dca,area1,area2
	# 0     1       2  3  4    5    6        7        8       9      10 11 12  13  14    15
	with open(infile) as fin:
		for line in fin:
			line=line.strip()
			if len(line) == 0:
				continue
			sarr = line.split(' ')
			pdb = sarr[0]
			chain = sarr[1]
			resi1 = sarr[2]
			resi2 = sarr[3]
			resn1 = sarr[4]
			resn2 = sarr[5]
			d1 = float(sarr[6])
			d2 = float(sarr[7])
			d3 = float(sarr[8])
			pfam = sarr[9]
			pos1 = sarr[10]
			pos2 = sarr[11]
			mip = float(sarr[12])
			dca = float(sarr[13])
			area1 = float(sarr[14])
			area2 = float(sarr[15])

			if (
					resn1 in cp.chargedaa and resn2 in cp.chargedaa and
					( d1 < 4.5 or d2 <4.5 or d3<4.5 ) and 
					( (area1 + area2) < 5 )
				):
				outlist.append(line)
	n = len(outlist)
	cp._info('%s %d records found.' % (infile, n))
	if n > 0:
		with open(outfile, 'w') as fout:
			fout.write('%s\n' % ('\n'.join(outlist)))

if __name__ == '__main__':
	cp.dispatch(__name__)