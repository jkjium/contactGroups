import sys
import numpy as np

# side-by-side rmsd stats
# input: .rmsd
#		cath174.pool.needle_EBLOSUM35_4_4.5.align.flat.rmsd
def main():
	if len(sys.argv) < 4:
		print 'Usage: python proc_rmsdstat.py cath174.pool.needle scsc2.3 b62'
		print 'output: cath174.pool.needle.scsc2.3_b62.rmsd.stat'
		return

	pref = sys.argv[1]
	matA = sys.argv[2]
	matB = sys.argv[3]

	rmsdfullA = []
	rmsdfullB = []
	print 'loading rmsd data ...'
	for i in xrange(2,18,2):
		for j in xrange(0, 9, 2):
			#print '%d-%d' % (i, j)
			rmsdfileA =  '%s_%s_%d_%.1f.align.flat.rmsd' % (pref, matA, i, j+0.5)
			rmsdfileB =  '%s_%s_%d_%.1f.align.flat.rmsd' % (pref, matB, i, j+0.5)

			with open(rmsdfileA) as pf:
				rmsdA = pf.readlines()

			with open(rmsdfileB) as pf:
				rmsdB = pf.readlines()

			if len(rmsdA)!=len(rmsdB):
				print 'error: unmatched flatfile'
				return

			rmsd_i_j_A = []
			rmsd_i_j_B = []
			for k in xrange(0, len(rmsdA)):
				arrA = rmsdA[k].split(' ')
				arrB = rmsdB[k].split(' ')

				nA = int(arrA[1])
				rA = float(arrA[2])

				nB = int(arrB[1])
				rB = float(arrB[2])

				# skip 0 rmsd entry
				if nA == 0 or nB == 0:
					continue

				rmsd_i_j_A.append(rA/nA)
				rmsd_i_j_B.append(rB/nB)

			rmsdfullA.append(rmsd_i_j_A)
			rmsdfullB.append(rmsd_i_j_B)

	print (len(rmsdfullA), len(rmsdfullB))
	# stats for 40x2 sets of rmsd
	outfile = '%s.%s_%s.rmsd.stat' % (pref, matA, matB)
	print 'save to %s' % outfile
	fout = open(outfile, 'w')
	for p in xrange(0, len(rmsdfullA)):
		totalA = np.sum(rmsdfullA[p])
		stdA = np.std(rmsdfullA[p])

		totalB = np.sum(rmsdfullB[p])
		stdB = np.std(rmsdfullB[p])

		fout.write('%f %f %f %f\n' % (totalA, stdA, totalB, stdB))
	fout.close()

if __name__ == '__main__':
	main()