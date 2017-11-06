import sys
import commp as cp

'''
#  Matrix made by matblas from blosum62.iij
#  * column uses minimum score
#  BLOSUM Clustered Scoring Matrix in 1/2 Bit Units
#  Blocks Database = /data/blocks_5.0/blocks.dat
#  Cluster Percentage: >= 62
#  Entropy =   0.6979, Expected =  -0.5209
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4
R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4
N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4
D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4
C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4
Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4
E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4
H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4
I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4
L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4
K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4
M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4
F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4
S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4
Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4
V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4
B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4
Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4
* -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1
'''

class smatrix(object):
	def __init__(self, smfile):
		scorelist = []
		with open(smfile) as fp:
			lines = fp.readlines()

		for line in lines:
			line = line.strip()
			if len(line) < 20: # at least 20 AA
				continue
			if line[0] == '#':
				continue
			if any(c.isdigit() for c in line):
				scorelist.append(line.split())
			else: # alphabet
				aa = line.split()
		#print aa
		# square matrix 'A:E' = 'E:A'
		self.score=dict( ('%s:%s' % (aa[i], aa[j-1]), int(scorelist[i][j])) for i in xrange(0,len(scorelist)) for j in xrange(1, len(aa)))
		self.aa = aa
		#print repr(self.score)

	def getscore(self, a1, a2):
		if a1 not in self.aa:
			cp._err('Invalid AA: %s' % a1)
		if a2 not in self.aa:
			cp._err('Invalid AA: %s' % a2)

		return self.score['%s:%s' % (a1,a2)]


	def psubdict(self):
		# first get the minimum and translate to remove all the negative values
		minscore = min(self.score.values())
		cp._info('min: %d' % minscore)
		order_aa = sorted([aa for aa in self.aa if aa not in cp.abaa])
		cp._info(order_aa)
		trans_sm = dict((k, self.score[k]-minscore) for k in self.score)
		cp._info(trans_sm)
		#for i in xrange(0,len(order_aa)):
		#	for j in xrange(i, len(order_aa)):




# for Dr. Jernigan 2017 Sep Grant
# mutant entry file as:
#	 Q04771	207	Q	E	
# $ python utils_sm.py batchscore b62 total.del.list mu.del.b62
def proc_batchscore():
	if len(sys.argv) < 4:
		cp._err('Usage: python utils_sm.py batchscore sm mutant_entyfile')
	sm = smatrix(sys.argv[2])

	with open(sys.argv[3]) as fp:
		mutant = [line.strip().split() for line in fp.readlines() if len(line) > 1]

	scores = [sm.getscore(m[2], m[3]) for m in mutant]
	#print repr(scores)

	'''
	outfile = sys.argv[4]
	with open(outfile, 'w') as fp:
		fp.write(', '.join([str(s) for s in scores]))
	print 'save to %s' % outfile
	'''
	print '= np.array([%s])' % ', '.join([str(s) for s in scores])


	'''
	for i in xrange(0, len(mutant)):
		print '%s %s %s %s : %2i' % (mutant[i][0], mutant[i][1], mutant[i][2], mutant[i][3], scores[i])
	'''


def psub():
	if len(sys.argv) < 2:
		cp._err('Usage: python utils_sm.py pusb b62 outfile')

	outfile = sys.argv[1]
	sm = smatrix(sys.argv[2])

	sm.psubdict()
	'''
	with open(outfile, 'w') as fp:
		sm.psubdict()
	'''


def test():
	'''
	sm = smatrix('seed7990_full.sm')
	print '%s:%s %d' % (sys.argv[2], sys.argv[3], sm.getscore(sys.argv[2], sys.argv[3]))
	'''
	pass


# main routine
def main():
	if len(sys.argv)<2:
		print 'Usage: python utils_sm.py cmd [args ...]'
		return

	dispatch = {
		'test':test,
		'batchscore': proc_batchscore,
		'psub': psub
	}

	cmd = sys.argv[1]

	if sys.argv[1] not in dispatch:
		print 'invalid cmd: %s' % sys.argv[1]
		return
	else:
		dispatch[sys.argv[1]]()

if __name__ == '__main__':
	main()
