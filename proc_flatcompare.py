import sys
import numpy as np

# ./needle 1hg8.seq 1mty.seq -gapopen 10 -gapextend 0.5 -data combined.matrix.0.002.new -aformat markx3 t.align
# parse emboss needle markx3 format output
# output a flat csv file
'''
$1:  file name
$2:  aligned sequence length
$3:  identity number
$4:  identity percentile 
$5:  similarity number
$6:  similarity percentile
$7:  gaps number
$8:  gaps percentile
$9:  align score
$10:  seq A pure length
$11: aligned seq A
$12: seq B pure length
$13: aligned seq B
'''


def main():
	if len(sys.argv) < 2:
		print 'python proc_flatcompare.py all-0.1.flat all-0.2..flat'
		print 'output: seq1-seq2.align.flat'
		return

	flat1 = sys.argv[1]
	flat2 = sys.argv[2]
	print '%s vs %s' % (flat1, flat2)

	f1 = open(flat1, 'r')
	lines1 = f1.readlines()
	f1.close()

	f2 = open(flat2, 'r')
	lines2 = f2.readlines()
	f2.close()

	identity = []
	similarity = []
	gaps = []
	for i in xrange(0, len(lines1)):
		line1 = lines1[i].strip()
		a1 = line1.split(' ')
		line2 = lines2[i].strip()
		a2 = line2.split(' ')
		if a1[0][0:17] != a2[0][0:17]:
			print 'error: unmatched titles: %s - %s' % (a1[0], a2[0])
			return
		identity.append(float(a1[3])-float(a2[3]))
		similarity.append(float(a1[5])-float(a2[5]))
		gaps.append(float(a1[7])-float(a2[7]))

	did=np.array(identity)
	dsm=np.array(similarity)
	dgp=np.array(gaps)

	print 'identity   +: %d\t-: %d\t.: %d' % ((did>0).sum(), (did<0).sum(), (did==0).sum())
	print 'similarity +: %d\t-: %d\t.: %d' % ((dsm>0).sum(), (dsm<0).sum(), (dsm==0).sum())
	print 'gaps       +: %d\t-: %d\t.: %d' % ((dgp>0).sum(), (dgp<0).sum(), (dgp==0).sum())




if __name__ == '__main__':
	main()
