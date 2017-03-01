import sys

# ./needle 1hg8.seq 1mty.seq -gapopen 10 -gapextend 0.5 -data combined.matrix.0.002.new -aformat markx3 t.align
# parse emboss needle markx3 format output
# output a flat csv file
'''
########################################
# Program: needle
# Rundate: Thu  1 Dec 2016 23:24:12
# Commandline: needle
#    [-asequence] 1hg8.seq
#    [-bsequence] 1mty.seq
#    -gapopen 10
#    -gapextend 0.5
#    -datafile combined.matrix.0.002.new
#    -aformat markx3
#    [-outfile] t.align
# Align_format: markx3
# Report_file: t.align
########################################

#=======================================
#
# Aligned_sequences: 2
# 1:
# 2:
# Matrix: combined.matrix.0.002.new
# Gap_penalty: 10.0
# Extend_penalty: 0.5
#
# Length: 451
# Identity:      76/451 (16.9%)
# Similarity:   125/451 (27.7%)
# Gaps:         186/451 (41.2%)
# Score: 75.5
#
#
#=======================================

> ..
------------------------CKNIVLNGFQVPTGKQL-DLSSLQ--
--------------NDSTVTFKG-------TTTFATTAD----ND-----
S
> ..
ERRRGLTDPEMAAVILKALPEAPLDGNNKMGYFVTPRWKRLTEYEALTVY
AQPNADWIAGGLDWGDWTQKFHGGRPSWGNETTELRTVDWFKHRDPLRRW
-

#---------------------------------------
#---------------------------------------
'''

def parseFasta(lines, i):
	AAset = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X', '-',
			 'a', 'r', 'n', 'd', 'c', 'q', 'e', 'g', 'h', 'i', 'l', 'k', 'm', 'f', 'p', 's', 't', 'w', 'y', 'v', 'b', 'z', 'x']
	fastalines = []
	while lines[i][0] in AAset:
		fastalines.append(lines[i].strip())
		i+=1
	return (''.join(fastalines), i-1)


def main():
	if len(sys.argv) < 2:
		print 'proc_needleParse: Too few parameters'
		print 'python proc_needleParse.py seq1-seq2.align'
		print 'output: seq1-seq2.align.flat'
		return

	alignfile = sys.argv[1]
	outfile = alignfile + '.flat'
	print 'alignment file: %s' % alignfile

	fin = open(alignfile, 'r')
	lines = fin.readlines()
	fin.close()

	i = 0
	out = []
	out.append(alignfile)

	while i<len(lines):
		#print lines[i]
		line = lines[i].strip()
		i+=1
		if len(line) < 1:
			continue

		# Gap_penalty: 10.0
		if 'Gap_penalty:' in line:
			 tsvline = ' '.join(line.split()) 
			 strArray = tsvline.split(' ')
			 gapopen = strArray[2]
			 out.append(gapopen)
			 print 'Gap_penalty: %s' % gapopen
			 continue

		# Extend_penalty: 0.5			 
		if 'Extend_penalty:' in line:
			 tsvline = ' '.join(line.split()) 
			 strArray = tsvline.split(' ')
			 gapextend = strArray[2]
			 out.append(gapextend)
			 print 'Extend_penalty: %s' % gapextend
			 continue			

		# Length: 451
		if 'Length:' in line:
			 tsvline = ' '.join(line.split()) 
			 strArray = tsvline.split(' ')
			 length = strArray[2]
			 out.append(length)
			 print 'length: %s' % length
			 continue

		# Identity:      76/451 (16.9%)	
		if 'Identity:' in line:
			 tsvline = ' '.join(line.split()) 
			 strArray = tsvline.split(' ')
			 ratioStr = strArray[2]
			 ratioArr = ratioStr.split('/')
			 out.append(ratioArr[0]) # number of identity
			 ratio = float(ratioArr[0])/float(ratioArr[1])
			 identity = '%.1f' % (ratio * 100)
			 out.append(identity)
			 print 'identity: %s' % identity
			 continue

		# Similarity:   125/451 (27.7%)
		if 'Similarity:' in line:
			 tsvline = ' '.join(line.split()) 
			 strArray = tsvline.split(' ')
			 ratioStr = strArray[2]
			 ratioArr = ratioStr.split('/')
			 out.append(ratioArr[0]) # number of similarity
			 ratio = float(ratioArr[0])/float(ratioArr[1])
			 similarity = '%.1f' % (ratio * 100)
			 out.append(similarity)
			 print 'similarity: %s' % similarity
			 continue

		# Gaps:         186/451 (41.2%)
		if 'Gaps:' in line:
			 tsvline = ' '.join(line.split()) 
			 strArray = tsvline.split(' ')
			 ratioStr = strArray[2]
			 ratioArr = ratioStr.split('/')
			 out.append(ratioArr[0]) # number of gaps
			 ratio = float(ratioArr[0])/float(ratioArr[1])
			 gaps = '%.1f' % (ratio * 100)
			 out.append(gaps)
			 print 'gaps: %s' % gaps
			 continue

		# Score: 75.5
		if 'Score:' in line:
			 tsvline = ' '.join(line.split()) 
			 strArray = tsvline.split(' ')
			 score = strArray[2]
			 out.append(score)
			 print 'score: %s' % score
			 continue

		# > ..
		if '>' == line[0]:
			fastaline, index = parseFasta(lines, i) # parse fasta sequence from the ith line
			out.append(('%d' % len(fastaline.replace('-', '')))) # pure sequence length
			out.append(fastaline)
			i = index 
			print 'seq: %s ... len: %d' % (fastaline[0:5], len(fastaline))
			continue

		#print repr(out)
	print 'save file: %s' % outfile
	fout = open(outfile, 'w')
	fout.write(' '.join(out)+'\n')
	fout.close()

if __name__ == '__main__':
	main()
