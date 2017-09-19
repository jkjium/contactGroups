import sys
import math
import commp as cp

def main():

	mfile = sys.argv[1]

	# each line.split() = ['Q04771', '375', 'R', 'P']
	mulist = set()
	with open(mfile) as fp:
		for line in fp:
			muarr = line.split()
			mulist.add((muarr[2], muarr[3]))

	#for m in mulist:
	#	print m

	print [math.fabs(cp.aaprop[m[0]][21] - cp.aaprop[m[1]][21]) for m in mulist]


if __name__ == '__main__':
	main()