import sys
def main():
	if len(sys.argv) < 3:
		print 'python gen_batchPairSeq.py seqfilelist matrix_name'
		print 'sequence file must in *.seq naming format'
		return

	fin = open(sys.argv[1])
	lines = fin.readlines()
	fin.close()

	matrix = sys.argv[2]

	out = []
	for i in xrange(0, len(lines)):
		for j in xrange(i+1, len(lines)):
			seq1 = lines[i].strip()
			seq2 = lines[j].strip()
			outfile = '%s-%s.%s.align' % (seq1[:-4], seq2[:-4], matrix)
			cmdstr = './needle %s %s -gapopen 10 -gapextend 0.5 -data %s -aformat markx3 %s' % (seq1, seq2, matrix, outfile)
			out.append(cmdstr)

	print '\n'.join(out)
	pass
if __name__ == '__main__':
	main()