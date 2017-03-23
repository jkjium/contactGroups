# to filter out the pairs that have similar correspondent in the pool
import sys
import subprocess as sp
from sets import Set

from alignflat import palign

def smcheck(s1a, s1b, s2a, s2b):
	#ret = sp.Popen(['aln.sh',s1a,s2a], stdout=sp.PIPE).communicate()[0]
	ret = sp.check_output(['aln.sh',s1a,s2a])
	pid1= float(ret[ret.index('(')+1:ret.index('%')])

	#ret = sp.Popen(['aln.sh',s1b,s2b], stdout=sp.PIPE).communicate()[0]
	ret = sp.check_output(['aln.sh',s1b,s2b])
	pid2= float(ret[ret.index('(')+1:ret.index('%')])

	if (pid1 >=90) and (pid2 >= 90):
		return False

	#ret = sp.Popen(['aln.sh',s1a,s2b], stdout=sp.PIPE).communicate()[0]
	ret = sp.check_output(['aln.sh',s1a,s2b])
	pid3= float(ret[ret.index('(')+1:ret.index('%')])

	#ret = sp.Popen(['aln.sh',s2a,s1b], stdout=sp.PIPE).communicate()[0]
	ret = sp.check_output(['aln.sh',s2a,s1b])
	pid4= float(ret[ret.index('(')+1:ret.index('%')])

	if (pid3 >=90) and (pid4 >= 90):
		return False

	# add into outlist
	return True

def main():
	if len(sys.argv)<3:
		print 'Usage: python proc_filterpair.py cathpair.flat outfile'
		exit(1)

	flatfile = sys.argv[1]
	outfile = sys.argv[2]

	with open(flatfile) as fp:
		lines = fp.readlines()

	outlist = []
	pa0 = palign(lines[0].strip())
	outlist.append(pa0)

	gap = ['.', '-']
	print len(lines)
	for i in xrange(1, len(lines)):
		flag = True
		#for j in xrange(0, len(outlist)):
		for p in outlist:
			pa = palign(lines[i].strip())
			s1a = p.seqA.translate(None, ''.join(gap))
			s1b = p.seqB.translate(None, ''.join(gap))
			s2a = pa.seqA.translate(None, ''.join(gap))
			s2b = pa.seqB.translate(None, ''.join(gap))

			#print '%s - %s' % (p.name, pa.name)
			if smcheck(s1a,s1b,s2a,s2b) == False:
				print 'rm: %s for %s' % (pa.name, p.name)
				flag = False
				break

		# no flase flag
		if flag == True:
			outlist.append(pa)

	# write outlist into file
	fout = open(outfile, 'w')			
	for p in outlist:
		fout.write(p.flatstr+'\n')
	fout.close()

	print '%d unique pairs saved to %s' % (len(outlist), outfile)

if __name__ == '__main__':
	main()