import sys
import numpy as np

def main():
	if len(sys.argv) < 3:
		print 'Usage: python proc_alignscore.py smfile flatfile'
		return
	
	smfile = sys.argv[1]
	flatfile = sys.argv[2]

	ab = {'A':0, 'R':1, 'N':2, 'D':3, 'C':4, 'Q':5, 'E':6, 'G':7, 'H':8, 'I':9, 'L':10, 'K':11, 'M':12, 'F':13, 'P':14, 'S':15, 'T':16, 'W':17, 'Y':18, 'V':19, 'B':20, 'Z':21, 'X':22, '*':23}

	# load matrix
	sm = np.loadtxt(smfile)

	# load flatfile	
	with open(flatfile) as f:
		flatlines = f.readlines()

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
	fout = open(flatfile+'.alignscore', 'w')
	for line in flatlines:
		line = line.strip()
		flatArray = line.split(' ')
		a1 = flatArray[10].upper()
		a2 = flatArray[12].upper()
		print flatArray[0] 
		#break
		for i in xrange(0, len(a1)):
			if a1[i]=='-' or a2[i] == '-': 
				continue
			if a1[i]!=a2[i]:
				print '%s - %s : %d' % (a1[i], a2[i], sm[ab[a1[i]]][ab[a2[i]]])
				fout.write('%s %s %s %d\n' % (flatArray[0], a1[i], a2[i], sm[ab[a1[i]]][ab[a2[i]]]))

	fout.close()



if __name__ == '__main__':
	main()