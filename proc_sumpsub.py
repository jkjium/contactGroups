import sys
import commp as cp
import collections

# q: quadstr
def quaddhat(q):
	return abs(int(cp.aaprop[q[0]][21])+int(cp.aaprop[q[1]][21]) - int(cp.aaprop[q[2]][21]) - int(cp.aaprop[q[3]][21]))


# iterate .allscol file to 
# sum up the pair substitution
def main():
	if len(sys.argv) < 2:
		cp._err('Usage: python proc_sumpsub.py 4-scol.allpsub')

	allpsubfile = sys.argv[1]
	psdict = collections.defaultdict(int)
	total = 0
	cp._info('loading psub ...')
	with open(allpsubfile) as fp:
		# ERNY 15096 t2
		for line in fp:
			sarr = line.strip().split(' ')
			psdict[sarr[0]]+=int(sarr[1])
			total+=float(sarr[1])

	outpsfile = allpsubfile + '.ps'
	cp._info('writing %s' % outpsfile)
	with open(outpsfile, 'w') as fp:
		for k in psdict:
			fp.write('%s %s %.6f %d\n' % (cp.quadtype(k), k, psdict[k]/total, quaddhat(k)))

	cp._info('save pair substitution: %s' % outpsfile)


if __name__ == '__main__':
	main()