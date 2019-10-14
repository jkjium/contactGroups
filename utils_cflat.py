'''
20190329 for finding charged buired amino acid pairs (coevolved)
[kjia@lhb-ps1 ~/workspace/pfam31.0/p90_ext] *.rcflat
'''

import commp as cp
import numpy as np

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
					d3 > 3 and d2 <4.5 and 
					( (area1 + area2) < 5 )
				):
				outlist.append(line)
	n = len(outlist)
	cp._info('%s %d records found.' % (infile, n))
	if n > 0:
		with open(outfile, 'w') as fout:
			fout.write('%s\n' % ('\n'.join(outlist)))

'''
append a ce value by the stub extracted from cflat
input: 
1. stub from cflat 
	msai1 msai2
	...
2. new ce value with msai pairs
	msai1 msai2 ce.value
	...
3. key columns
	indicate the number of columns are used as a key to match the row from stubfile
	columns for key must be columns at front of a cefile
	rest of the columns are output together
output:
	values from .cefile corresponding to the stub order
	print msai for debug
'''
def appendbystub(arglist):
	if len(arglist) < 4:
		cp._err('Usage: utils_cflat.py appendbystub cflat.stub.file new.ce.file 2 out.mapped.ce.vec')

	stubfile = arglist[0]
	cefile = arglist[1]
	knum = int(arglist[2])
	outfile = arglist[3]

	# load ce into a dictionary
	cedict = {}
	celist = cp.loadlines(cefile)
	cp._info('%d ce values loaded' % len(celist))
	for line in celist:
		sarr = line.split(' ')
		key = ' '.join(sarr[:knum])
		value = ' '.join(sarr[knum:])
		cedict[key] = value

	fill = ' '.join(['-191']*(len(line.split(' '))-knum))
	# map stub from cedict
	outlist = []
	count = 0
	for msaipair in cp.loadlines(stubfile):
		if msaipair in cedict:
			outlist.append(cedict[msaipair])
		else:
			cp._info("key %s has no match in celist" % msaipair)
			count+=1
			outlist.append(fill)
	cp._info('%d records mapped' % len(outlist))
	cp._info('%d records missed' % count)

	# save outfile
	with open(outfile, 'w') as fout:
		fout.write('%s\n' % '\n'.join(outlist))
	cp._info('save to %s' % outfile)

# calculate column based zscore for each elements
def zscore(arglist):
	if len(arglist) < 2:
		cp._err('Usage: utils_cflat.py in_value_matrix.file out_zscore_matrix.file')
	infile = arglist[0]
	outfile = arglist[1]
	d = np.loadtxt(infile, delimiter=' ')
	# d.shape = (row, column)
	zscorelist = [cp.zscore(d[:,i]) for i in xrange(0, d.shape[1])]
	np.savetxt(outfile, np.array(zscorelist).T, fmt='%.4f', delimiter=' ')
	cp._info('save zscore to %s' % outfile)







































if __name__ == '__main__':
	cp.dispatch(__name__)