import sys

def main():
	if len(sys.argv) < 3:
		print 'Usage: python proc_1sthead.py PF00000_full.txt ec|w'
		exit()

	msafile = sys.argv[1]
	opt = sys.argv[2]
	with open(msafile) as fp:
		for line in fp:
			strarr = line[1:].split('/')
			# calculate_evolutionary_constraints('PF00071_v25.fa','RASH_HUMAN','your_output_filename.txt');
			if opt == 'ec':
				print 'calculate_evolutionary_constraints(\'%s\',\'%s\',\'%s.dca\');' % (msafile, strarr[0], msafile)
			else:
				print 'evfold_weight(\'%s\',\'%s\', 0.3);' % (msafile, strarr[0])
			break

if __name__ == '__main__':
	main()