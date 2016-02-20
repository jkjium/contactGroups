# output contact group in WenZhou's format
import sys
import itertools
from cg import cg
from naccess import naccess

# varset: list of alphabet sets
# varset = [list1, list2, list3]
# [['X', 'Y', 'Z'], ['B', 'E'], ['1', '2']]
# output: ['XB1', 'XB2', 'XE1', 'XE2', 'YB1', 'YB2', 'YE1', 'YE2', 'ZB1', 'ZB2', 'ZE1', 'ZE2']
def expandVars(varset):
	return [''.join(ll) for ll in itertools.product(*varset)]


def main():
	if len(sys.argv) < 3:
		print "Usage python proc_scoreboard.py cg_file rsa_file"
		return

	AAlist = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
	cgfile = sys.argv[1]
	nafile = sys.argv[2]
	outfile = cgfile+'.nascore'

	print 'loading %s, %s' % (cgfile, nafile)

	na = naccess(nafile)
	alphabet = expandVars([AAlist, na.alphabet])
	#print repr(alphabet)

	cgs = [cg(line.strip(), alphabet) for line in open(cgfile)]

	fo = open(outfile, 'w')
	for c in cgs:
		if len(c.AAgroup) > 1:
			fo.write(c.nascore(na)+'\n')
	fo.close()

if __name__=="__main__":
	main()

