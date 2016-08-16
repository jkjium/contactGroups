# output contact group in contact group score + volume + rsa
import sys
import itertools
from cg import cg
from naccess import naccess

def main():
	if len(sys.argv) < 3:
		print "Usage python proc_nvboard.py cg_file rsa_file"
		print "python proc_nvboard.py 1k2p.tip.hcg 1k2p.rsa"
		print "output: 1k2p.tip.hcg.nvscore"
		return

	AAlist = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
	cgfile = sys.argv[1]
	nafile = sys.argv[2]
	outfile = cgfile+'.nvscore'

	#print 'loading %s, %s' % (cgfile, nafile)

	na = naccess(nafile)

	cgs = [cg(line.strip(), '') for line in open(cgfile)]

	fo = open(outfile, 'w')
	for c in cgs:
		if len(c.AAgroup) > 1:
			fo.write(c.nvscore(na)+'\n')
	fo.close()
	print 'finish writing %s.' % outfile

if __name__=="__main__":
	main()

