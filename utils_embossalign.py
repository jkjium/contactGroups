import sys
import subprocess

'''
new version of alignflat and proc_testcase.py

'''
class embossalign(object):
	def __init__(self, flatStr):
		flatArray = flatStr.split(' ')
		if len(flatArray)!= 17:
			print 'Invalid flat string len: %d\n' % len(flatArray)
			print flatStr
			exit(1)
		self.flatstr = flatStr
		self.name = flatArray[0]
		self.program = flatArray[1]
		self.matrix = flatArray[2]
		self.gapopen = float(flatArray[3])
		self.gapextend = float(flatArray[4])
		self.alignlen = int(flatArray[5])
		self.nid = float(flatArray[6])
		self.pid = float(flatArray[7])
		self.nsm = float(flatArray[8])
		self.psm = float(flatArray[9])
		self.ngp = float(flatArray[10])
		self.pgp = float(flatArray[11])
		self.score = float(flatArray[12])
		self.seqAlen = float(flatArray[13])
		self.seqA = flatArray[14]
		self.seqBlen = float(flatArray[15])
		self.seqB = flatArray[16]


	# dump object to stdout
	def dump(self):
		print '\n-----------------------------'
		print 'name: %s' % self.name
		print 'program: %s' % self.program
		print 'matrix: %s' % self.matrix
		print 'gap open penalty: %f' % self.gapopen
		print 'gap extend penalty: %f' % self.gapextend
		print 'alignment length: %d' % self.alignlen
		print 'number of identity: %d' % self.nid
		print 'percent of identity: %f' % self.pid
		print 'number of similarity: %d' % self.nsm
		print 'percent of similarity: %f' % self.psm
		print 'number of gaps: %d' % self.ngp
		print 'percent of gaps: %f' % self.pgp
		print 'alignment score: %f' % self.score
		print 'seq A length: %d' % self.seqAlen
		print 'aligned seq A: %s' % self.seqA
		print 'seq B length: %d' % self.seqBlen
		print 'aligned seq B: %s' % self.seqB
		print '-----------------------------\n'


# call needle or water to get the alignment 
# needle <(echo -e ">a\n$1") <(echo -e ">b\n$2") "${@:3}" -filter
def align_exec(s1, s2, cmd='needle', matrix='B62', gapopen='10.0', gapextend='0.5'):
	#$ ./align.sh needle "ALIGN" "LINE" B62 8 2
	if (('%d.%d') % (sys.version_info.major, sys.varsion_info.minor)) == '2.6':
		ret = subprocess.Popen(['align.sh', cmd, s1, s2, matrix, gapopen, gapextend], stdout=subprocess.PIPE).communicate()[0].strip()
	else:	
		ret = subprocess.check_output(['align.sh', cmd, s1, s2, matrix, gapopen, gapextend])
	return ret
 

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


def help(d):
	print '--------------------'
	for name in d:
		d[name]([])
		print ''
	print '--------------------'

def main():
	dispatch = {
		'printpair':printpair,
		'testpool':testpool
	}

	if len(sys.argv) < 2:
		help(dispatch)
		return
	cmd = sys.argv[1]
	if cmd not in dispatch:
		print 'Invalid cmd %s' % cmd
		return

	dispatch[sys.argv[1]](sys.argv)

if __name__ == '__main__':
	main()