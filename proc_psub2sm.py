import sys
import math
from collections import defaultdict
import commp as cp

'''
$ sort -r -gk 4 4-tip-scol.allpsub.norm.t2|head
t2 DREK 92688902 7436.230808 0.248543 4472289669945.905273 336046726.166748 1
t2 DKER 65872545 6530.121776 0.185251 3178385935961.856445 238823123.522618 3
t2 DSET 67216576 5687.592443 0.103695 1853011605699.154297 214848082.010815 2
'''
def main():
	if len(sys.argv) < 3:
		print 'Usage: python proc_psub2sm.py psubfile outfile'
		exit()

	psubfile = sys.argv[1]
	outfile = sys.argv[2]

	p = defaultdict(float)
	true_psub = {}
	total = 0.0
	with open(psubfile) as fp:
		for line in fp:
			line = line.strip()
			if len(line)==0:
				continue
			strarr = line.split(' ')
			f = float(strarr[3])
			k=strarr[1]
			true_psub[k] = f

			# DREK = 7436.230808
			# DE = ED = 7436.230808
			k1 = '%s%s' % (k[0],k[2])
			k12 = '%s%s' % (k[2], k[0])
			p[k1] += f
			p[k12] += f

			# RK = KR = 7436.230808
			k2 = '%s%s' % (k[1],k[3])
			k22 = '%s%s' % (k[3],k[1])
			p[k2] += f
			p[k22] += f

			total+=f

	# process true psub into log
	for tk in true_psub:
		true_psub[tk] = math.log(true_psub[tk]/total, 2)
	minscore = min(true_psub.values())
	for tk in true_psub:
		true_psub[tk]-=minscore

	# process pseudo psub into log 
	for k in p:
		p[k] = math.log(p[k]/total,2)
	trans_sm = sorted([(k, p[k]-min(p.values())) for k in p], key=lambda x: x[1], reverse=True)

	#trans_p = dict((k, v) for k,v in trans_sm)
	psepsub = {}
	for k,v in trans_sm:
		for k1,v1 in trans_sm:
			psubstr = cp.quad_permu([k[0],k1[0],k[1],k1[1]])
			t = cp.quadtype(psubstr)
			if t == 't2' or t == 't21':
				#print (psubstr, int(round(v)), int(round(v1)), int(round(v+v1)))
				psepsub[psubstr] = v+v1

	# extract all key from true psub
	norm_psepsub = cp.rank01(dict((k, psepsub[k]) for k in true_psub))
	norm_truepsub = cp.rank01(true_psub)

	with open(outfile, 'w') as fp:
		for k in true_psub:
			diff = norm_truepsub[k]-norm_psepsub[k]
			fp.write('%s %.8f %.8f %.8f %s\n' % (k, norm_truepsub[k], norm_psepsub[k], math.fabs(diff), ('-' if diff >=0 else '+')))

if __name__ == '__main__':
	main()	