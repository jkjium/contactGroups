import sys
import subprocess

'''
Matching a given sequence to a sequence pool
sequence pool in flat sequence format:
[length] [seqName] [sequence]
'''
def fetchflatseq(seqline):
	seqArray = seqline.split(' ')
	if len(seqArray)!=3:
		print 'invalid sequence line found.\n%s' % seqline
		return
	return (seqArray[1], seqArray[2].strip())


def main():
	if len(sys.argv) < 5:
		print 'Usage: python proc_seqmatch.py target seqpool {local, global} {B62, SU}'
		print 'Usage: python proc_seqmatch.py 4msp.flatseq seqpool.flat local B62'
		print 'target.seq is in flat format as well'
		exit(1)

	for p in sys.argv:
		print p,
	print ''

	targetfile = sys.argv[1]
	seqpoolfile = sys.argv[2]
	method = sys.argv[3]
	sm = sys.argv[4]

	cmd = './global.sh'

	if method == 'local':
		cmd = './local.sh'


	with open(targetfile) as fp:
		targetName, targetSeq = fetchflatseq(fp.readline())

	#print 'Taget name: %s' % targetName
	#print 'Taget seq:  %s' % targetSeq

	maxpid = 0.0
	with open(seqpoolfile) as fp:
		for line in fp:
			name, seq = fetchflatseq(line)
			#print '%s\n%s' % (name, seq)
			# output = subprocess.Popen(['ls', '-l'], stdout=subprocess.PIPE).communicate()[0]
			# ret = sp.check_output(['./aln.sh',si,sj])
			ret = subprocess.Popen([cmd, targetSeq, seq, '-data', sm], stdout=subprocess.PIPE).communicate()[0].strip()
			strpid= ret[ret.index('(')+1:ret.index('%')].strip()
			pid = float(strpid)			

			rettsv = ' '.join(ret.split())
			idstrArray = rettsv.split(' ')
			nid = int(idstrArray[2][0:idstrArray[2].index('/')])


			print '%s %d %.1f %d %.1f' % (name, len(seq), pid, nid, 100.0*nid/len(targetSeq))
			'''
			if int(strpid[0]) < 10:
				strpid = '0'
			outname = 'p.%s.%s.%s.%s' % (targetName,name,sm,strpid[0])
			fout = open(outname, 'w')
			fout.write('%s %s' % (targetSeq, seq))
			fout.close()		
			'''

if __name__ == '__main__':
	main()