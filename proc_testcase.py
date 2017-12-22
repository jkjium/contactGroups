import sys
import subprocess as sp 

from alignflat import palign


'''
def testpool():
	if len(sys.argv)<4:
		print 'Usage: python proc_testcase.py aligndiff testpool.list'
		return


	newsm = sys.argv[2]
	pairlist = sys.argv[3]

	pairarray = []
	with open(pairlist) as fp:
		for line in fp:
			pairarray.append(line.strip())
	print '%d pairs test cases loaded ..' % len(pairarray)


def aligndiff(pairpool, cmd='needle', gapopen='10.0', gapextend='0.5'):

	for p in pairpool:
		ret = align_exec(p, cmd, gapopen, gapextend)

'''

# parse FASTA sequence from markx3 align outputs
def parseFasta(lines, i):
	AAset = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X', '-',
			 'a', 'r', 'n', 'd', 'c', 'q', 'e', 'g', 'h', 'i', 'l', 'k', 'm', 'f', 'p', 's', 't', 'w', 'y', 'v', 'b', 'z', 'x']
	fastalines = []
	while lines[i][0] in AAset:
		fastalines.append(lines[i])
		i+=1
	return (''.join(fastalines), i-1)


# parse markx3 format into flat format
def alignparse(title, alignstr):
	lines = filter(None, alignstr.split('\n'))

	i = 0
	out = []
	out.append(title)

	while i<len(lines):
		#print lines[i]
		line = lines[i].strip()
		i+=1

		if 'Program:' in line:
			 tsvline = ' '.join(line.split()) 
			 strArray = tsvline.split(' ')
			 program = strArray[2]
			 out.append(program)
			 #print 'Program: %s' % program
			 continue

		# Matrix: SU
		if 'Matrix:' in line:
			 tsvline = ' '.join(line.split()) 
			 strArray = tsvline.split(' ')
			 matrix = strArray[2]
			 out.append(matrix)
			 #print 'Matrix: %s' % matrix
			 continue

		# Gap_penalty: 10.0
		if 'Gap_penalty:' in line:
			 tsvline = ' '.join(line.split()) 
			 strArray = tsvline.split(' ')
			 gapopen = strArray[2]
			 out.append(gapopen)
			 #print 'Gap_penalty: %s' % gapopen
			 continue

		# Extend_penalty: 0.5			 
		if 'Extend_penalty:' in line:
			 tsvline = ' '.join(line.split()) 
			 strArray = tsvline.split(' ')
			 gapextend = strArray[2]
			 out.append(gapextend)
			 #print 'Extend_penalty: %s' % gapextend
			 continue			

		# Length: 451
		if 'Length:' in line:
			 tsvline = ' '.join(line.split()) 
			 strArray = tsvline.split(' ')
			 length = strArray[2]
			 out.append(length)
			 #print 'length: %s' % length
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
			 #print 'identity: %s' % identity
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
			 #print 'similarity: %s' % similarity
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
			 #print 'gaps: %s' % gaps
			 continue

		# Score: 75.5
		if 'Score:' in line:
			 tsvline = ' '.join(line.split()) 
			 strArray = tsvline.split(' ')
			 score = strArray[2]
			 out.append(score)
			 #print 'score: %s' % score
			 continue

		# > ..
		if '>' == line[0]:
			fastaline, index = parseFasta(lines, i) # parse fasta sequence from the ith line
			out.append(('%d' % len(fastaline.replace('-', '')))) # pure sequence length
			out.append(fastaline)
			i = index 
			#print 'seq: %s ... len: %d' % (fastaline[0:5], len(fastaline))
			continue

	return ' '.join(out)


# call needle or water to get the alignment 
# needle <(echo -e ">a\n$1") <(echo -e ">b\n$2") "${@:3}" -filter
def align_exec(seqpair, cmd='needle', matrix='B62', gapopen='10.0', gapextend='0.5'):
	with open(seqpair) as fp:
		line = fp.readline()
		if len(line)<3:
			print 'Error: invalid pairfile : %s\n%s' % (seqpair, line)
			exit(1)
		seqArr = line.split(' ')

	#$ ./align.sh needle "ALIGN" "LINE" B62 8 2
	return sp.Popen(['align.sh', cmd, seqArr[0], seqArr[1], matrix, gapopen, gapextend], stdout=sp.PIPE).communicate()[0]
	#return sp.check_output(['align.sh', cmd, seqArr[0], seqArr[1], matrix, gapopen, gapextend])
 

def printpair():
	if len(sys.argv)< 3:
		print 'Usage: python proc_testcase.py alignpair pairfile B62'
		return
    
	pairfile = sys.argv[2]
	matrix = 'EBLOSUM62'
	if len(sys.argv) == 4:
		matrix = sys.argv[3]

	ret = align_exec(pairfile, 'needle', matrix, '10.0', '0.5')
	ap = palign(alignparse(pairfile, ret))
	ap.dump()


# output all flat files 
def alignpool(seqpool, param):
	#cathpair88.pool.needle_B62_8_8.5.align.flat
	flatfile = '%s.%s_%s_%s-%s_.flat' % (seqpool, param[0], param[1], param[2], param[3])
	#print 'save to: %s' % flatfile
	identity = 0
	similarity = 0
	gaps = 0


	fout = open(flatfile, 'w')
	with open(seqpool) as fp:
		for line in fp:
			name = line.strip()
			ret = align_exec(name, param[0], param[1], param[2], param[3]) # needle, B62, 10, 0.5
			flat = alignparse(name, ret)
			fout.write('%s\n' % flat)

			ap = palign(flat)
			identity+=ap.nid
			similarity+=ap.nsm
			gaps+=ap.ngp

	fout.close()
	print '%d %d %d %s %s %s %s %s' % (identity, similarity, gaps, seqpool, param[0], param[1], param[2], param[3])


# output all flat files 
def alignpoolshow(seqpool, param):
	#cathpair88.pool.needle_B62_8_8.5.align.flat
	flatfile = '%s.%s_%s_%s-%s_.flat' % (seqpool, param[0], param[1], param[2], param[3])
	#print 'save to: %s' % flatfile
	identity = 0
	similarity = 0
	gaps = 0


	#fout = open(flatfile, 'w')
	with open(seqpool) as fp:
		for line in fp:
			name = line.strip()
			ret = align_exec(name, param[0], param[1], param[2], param[3]) # needle, B62, 10, 0.5
			flat = alignparse(name, ret)
			#fout.write('%s\n' % flat)

			ap = palign(flat)
			identity+=ap.nid
			similarity+=ap.nsm
			gaps+=ap.ngp

	#fout.close()
	print '%d %d %d %s %s %s %s %s' % (identity, similarity, gaps, seqpool, param[0], param[1], param[2], param[3])


def testpool():
	if len(sys.argv)< 6:
		print 'Usage: python proc_testcase.py foo pairlist.txt needle B62 10.0 0.5'
		exit(1)

	pairlist = sys.argv[2]
	cmd = sys.argv[3]
	matrix = sys.argv[4]
	gapopen = sys.argv[5]
	gapextend = sys.argv[6]

	alignpool(pairlist, (cmd, matrix, gapopen, gapextend))


# same as testpool without writing flat file
def showpool():
	if len(sys.argv)< 6:
		print 'Usage: python proc_testcase.py foo pairlist.txt needle B62 10.0 0.5'
		exit(1)

	pairlist = sys.argv[2]
	cmd = sys.argv[3]
	matrix = sys.argv[4]
	gapopen = sys.argv[5]
	gapextend = sys.argv[6]

	alignpoolshow(pairlist, (cmd, matrix, gapopen, gapextend))



def main():
	dispatch = {
		#	'testpool':testpool, 
		'printpair':printpair,
		'testpool':testpool

	}

	if len(sys.argv)<3:
		for k in dispatch:
			dispatch[k]()
		exit(1)

	#print repr(sys.argv)

	cmd = sys.argv[1]

	flag = False
	for key in dispatch:
		if key == cmd:
			dispatch[key]()
			flag = True
	if flag == False:
		print 'Wrong cmd string'

if __name__ == '__main__':
	main()