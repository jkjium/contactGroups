# output contact group in WenZhou's format
import sys
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

	print 'loading %s, %s' % (cgfile, nafile)

	na = naccess(nafile)
	alphabet = expandVars([AAlist, na.alphabet])

	cgs = [cg(line.strip(), alphabet) for line in open(cgfile)]

	for c in cgs:
		

	fin = open(, 'r')
	fo = open(sys.argv[1]+'.nascore', 'w')
	for line in fin.readlines():
		cg = cgroup(line.strip())
		fo.write(cg.scoreboard2str()+'\n')
	fin.close()
	fo.close()

	'''
	#arrSingleVar = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
	#arrDoubleVar = ['AA','CC','DD','EE','FF','GG','HH','II','KK','LL','MM','NN','PP','QQ','RR','SS','TT','VV','WW','YY']
	arrTripleVar = ['AAA','CCC','DDD','EEE','FFF','GGG','HHH','III','KKK','LLL','MMM','NNN','PPP','QQQ','RRR','SSS','TTT','VVV','WWW','YYY']

	fin = open(sys.argv[1], 'r')
	out_name = sys.argv[1]+'.score'
	fout = open(out_name, 'w')
	for line in fin.readlines():
		line = line.strip()
		strArr = line.split(',')
		if len(strArr[1])==1:
			continue
		groupStr = ''.join(sorted(strArr[1]))
		board = ''
		for i in xrange(0,len(arrTripleVar)):
			if arrTripleVar[i] in groupStr:
				outStr = '0 0 1'
			elif arrTripleVar[i][0:2] in groupStr:
				outStr = '0 1 0'
			elif arrTripleVar[i][0] in groupStr:
				outStr = '1 0 0'
			else:
				outStr = '0 0 0'
			board = ('%s %s') % (board, outStr)

		board = board.lstrip(' ')
		fout.write(('%s,%s,%s,%s\n') % (strArr[0], board, strArr[1], strArr[2]))

	fin.close()
	fout.close()
	'''
if __name__=="__main__":
	main()

