import sys
import numpy as np

# side-by-side rmsd stats
def main():
	if len(sys.argv) < 3:
		print 'Usage: python proc_rmsdstat.py cath174.b62 cath174.cb2'
		return

	prefA = sys.argv[1]
	prefB = sys.argv[2]

	rmsdfullA = []
	rmsdfullB = []
	print 'loading rmsd data ...'
	for i in xrange(2,18,2):
		for j in xrange(0, 9, 2):
			print '%d-%d' % (i, j)
			rmsdfileA =  '%s.%d-%.1f.flat.rmsd' % (prefA, i, j+0.5)
			rmsdfileB =  '%s.%d-%.1f.flat.rmsd' % (prefB, i, j+0.5)

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
	print 'save to %s, %s' % (prefA+'.stat', prefB+'.stat')
	foutA = open(prefA+'.rmsd.stat', 'w')
	foutB = open(prefB+'.rmsd.stat', 'w')
	for p in xrange(0, len(rmsdfullA)):
		totalA = np.sum(rmsdfullA[p])
		stdA = np.std(rmsdfullA[p])
		foutA.write('%f %f\n' % (totalA, stdA))

		totalB = np.sum(rmsdfullB[p])
		stdB = np.std(rmsdfullB[p])
		foutB.write('%f %f\n' % (totalB, stdB))

	foutA.close()
	foutB.close()

if __name__ == '__main__':
	main()