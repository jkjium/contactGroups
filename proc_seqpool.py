import sys
import subprocess as sp 

def fetchseq(db, i):
	ci = '%dq;d' % i 
	ret = sp.check_output(['sed',ci, db]).strip()
	retArray = ret.split(' ')
	ni = retArray[1]
	si = retArray[2]	
	return (ni, si)


def main():
	if len(sys.argv)<5:
		print 'Usage: python proc_seqpool.py pfam-A-single.fasta 30-100 9347981 60000 0'
		print 'Usage: python proc_seqpool.py 2m.4m.flat 30-100 2000000 50000 1'
		return

	for p in sys.argv:
		print p,
	print ''

	dbname = sys.argv[1]
	cutoffArray = sys.argv[2].split('-')
	if len(cutoffArray)!=2:
		print 'invalid search range: %s' % sys.argv[2]
		exit()
	cmin = float(cutoffArray[0])
	cmax = float(cutoffArray[1])
	total = int(sys.argv[3])
	numlimit = int(sys.argv[4])
	bucket = int(sys.argv[5]) 


	print 'searching db:    %s' % dbname
	print 'searching range: %%%.1f - %%%.1f' % (cmin, cmax)
	print 'total seq:       %d' % total
	print 'pair limit:      %d' % numlimit

	seqpool = []
	if bucket == 1:
		print 'loading seqpool ...'
		sys.stdout.flush()
		with open(dbname) as fp:
			for line in fp: # line format: len seqName sequence
				strArray = line.split(' ')
				seqpool.append((strArray[1], strArray[2].strip()))
		print '%d sequences loaded.' % len(seqpool)


	# init counters ...
	count=0
	tag=1
	forwardflag = False

	if bucket == 1:
		i = 0
		j = 1
		ni = seqpool[i][0]
		si = seqpool[i][1]
	else:
		i = 1
		j = 2
		ni,si = fetchseq(dbname, i)


	# main loop
	while(True):
		if forwardflag == True:
			i+=1
			print 'forwarding i: %d' % i
			j=i+1
			if bucket == 1:
				ni = seqpool[i][0]
				si = seqpool[i][1]
			else:
				(ni,si) = fetchseq(dbname, i)

		if i == total:
			print 'Primary index i reaches the end. process halt.'
			exit()

		if j == total:
			forwardflag = True
			print 'j reach end. forward i.'
			continue

		if bucket == 1:
			nj = seqpool[j][0]
			sj = seqpool[j][1]
		else:
			nj,sj = fetchseq(dbname, j)


		#r = 1.0*abs(len(si)-len(sj))/min(len(si),len(sj))
		#print 'r: %.3f' % r

		#if r <= 0.3: 
		ret = sp.check_output(['./aln.sh',si,sj])
		#print ret
		strpid= ret[ret.index('(')+1:ret.index('%')]
		pid = float(strpid) # get identity percentage
		strpid = '%02d' % pid
		#print '[%d] - i: %d j: %d pid: %.2f ' % (count, i,j,pid)
		#print 'i: %d len(si): %d ni: %s\n%s' % (i, len(si), ni, si)
		#print 'j: %d len(sj): %d nj: %s\n%s\n' % (j, len(sj), nj, sj)

		# writing pair into file
		if pid <= cmax and pid > cmin :
			outname = 'p.%s.%s.%s' % (ni,nj,strpid[0])
			print '[%d] - i: %d j: %d pid: %.2f - [%s]\n' % (count, i,j,pid,outname)

			fout = open(outname, 'w')
			fout.write('%s %s' % (si, sj))
			fout.close()		

			count+=1
		else:
			# not in the identity range
			# keep i unchanged
			forwardflag = False
			#print '%d %d pid: %.1f' % (i, j, pid)
			#tag+=1
			#if tag%100==0:
			#	print '.',
			#	sys.stdout.flush()
			'''
			print 'i: %d len(si): %d ni: %s\n%s' % (i, len(si), ni, si)
			print 'j: %d len(sj): %d nj: %s\n%s\n' % (j, len(sj), nj, sj)
			'''
		'''
		else: # not in the same length level
			# forward i
			forwardflag = True
			tag+=1
			if tag%100==0:
				print '+',
				sys.stdout.flush()
			print '%d %d r: %.1f' % (i, j, r)
		'''

		j+=1
		if count == numlimit:
			print 'collection finished.'
			exit()


if __name__ == '__main__':
	main()
