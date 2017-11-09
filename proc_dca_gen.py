import sys

def main():
	if len(sys.argv) < 2:
		print 'Usage: python proc_1sthead.py PF00000_full.txt'
		exit()

	msafile = sys.argv[1]
	with open(msafile) as fp:
		for line in fp:
			strarr = line[1:].split('/')
			# calculate_evolutionary_constraints('PF00071_v25.fa','RASH_HUMAN','your_output_filename.txt');
			print 'fprintf \'%s\n\';calculate_evolutionary_constraints(\'%s\',\'%s\',\'%s.dca\');' % (msafile, msafile, strarr[0], msafile)
			# evfold_weight('PF00589_full.txt','A0A010NMT9_9MICC',0.3);
			print 'fprintf \'%s\n\';evfold_weight(\'%s\',\'%s\',0.3);' % (masfile, msafile, strarr[0])
			break

if __name__ == '__main__':
	main()