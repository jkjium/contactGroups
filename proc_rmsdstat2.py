import sys
import numpy as np
import commp as cp

# side-by-side rmsd stats for one .rmsd file
# [kjia@lhb-ps1 ~/workspace/pfam31.0/testcase/cath.testpool/selected.flatrmsd] 
'''
p.1a0pA02.1ae9A00.1-10-443-10 162 7.7363
p.1a0pA02.1aihA00.1-10-443-10 154 9.3331
p.1a0pA02.1xo0A02.1-10-443-10 162 15.5358
p.1a0pA02.2a3vA02.1-10-443-10 165 13.1172
'''
# for scsc3 re-submit
# print penalty, mean_rmsd, std
def main():
	if len(sys.argv) < 2:
		cp._err('Usage: python proc_rmsdstat_scsc3.py cath473.pool.needle_B62_6-3.5_.flat.rmsd')

	infile = sys.argv[1]
	sarr = infile.split('_')
	penalty = sarr[2]
	sarr2 = penalty.split('-')
	index = float(sarr2[0])*10+float(sarr2[1])

	mean_rmsd = [] 
	with open(infile) as fp:
		for line in fp:
			line = line.strip()
			if len(line) == 0:
				continue
			sarr = line.split(' ')

			nres = float(sarr[1])
			rmsd = float(sarr[2])
			if nres >0:
				mean_rmsd.append(rmsd/nres)
			else:
				mean_rmsd.append(0.0)

	nrmsd = np.array(mean_rmsd)
	print '%s %.4f %.4f %.2f' % (penalty, nrmsd.mean(), nrmsd.std(), index)

if __name__ == '__main__':
	main()