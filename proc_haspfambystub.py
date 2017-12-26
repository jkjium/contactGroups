import sys
import commp as cp
import collections
from utils_pfamscan import utils_pfamscan

def main():
	if len(sys.argv) < 4:
		cp._err('Usage: python proc_haspfambystub.py stubfile jsonfilelistfile tpfpoutfile')

	stubfile = sys.argv[1]
	jsonfilelistfile = sys.argv[2]
	outfile =sys.argv[3]

	# load stub file
	# uc0003 PF00000
	with open(stubfile) as fp:
		pfamstub=dict((line[0:6], line[7:14]) for line in fp if len(line.strip())!=0)
	cp._info('%d stub loaded.' % len(pfamstub))

 	tp = collections.defaultdict(int)
 	fp = collections.defaultdict(int)
	with open(jsonfilelistfile) as fp:
		for line in fp:
			line = line.strip()
			if len(line) == 0:
				continue
			# uc0005.nr.b62.uniq.out0000.fa.json
			sarr = line.split('.')
			pfamid = pfamstub[sarr[0]]
			ups = utils_pfamscan(line)
			if ups.getMatchpfs(pfamid):
				# uc0005.b62
				print '%s.%s %s 1' % (sarr[0], sarr[2], pfamid)
				tp['%s.%s' % (sarr[0], sarr[2])]+=1
			else:
				print '%s.%s %s 0' % (sarr[0], sarr[2], pfamid)
				fp['%s.%s' % (sarr[0], sarr[2])]+=1

	with open(outfile, 'w') as fp:
		for k in pfamstub:
			b62tp = tp['%s.b62' % k]
			b62fp = fp['%s.b62' % k]
			scsctp = tp['%s.scsc' % k]
			scscfp = fp['%s.scsc' % k]
			# uc0003 PF00001 b62.tp b62.fp scsc.tp scsc.fp
			fp.write('%s %s %d %d %d %d\n' % (k, pfamstub[k], b62tp, b62fp, scsctp, scscfp))


if __name__ == '__main__':
	main()