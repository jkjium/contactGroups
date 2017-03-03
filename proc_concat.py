import sys
from msa import msa

def main():
	if len(sys.argv)< 4:
		print 'Usage: python proc_concat.py Dsg2_aligned.fasta Ecad_aligned.fasta pairfile outfile'
		exit()

	msa1 = sys.argv[1]
	msa2 = sys.argv[2]
	pairfile = sys.argv[3]
	outfile = sys.argv[4]

	m1 = msa(msa1)
	m2 = msa(msa2)

	print '%s len: %d, num: %d' % (msa1, m1.seqlen, m1.seqNum)
	print '%s len: %d, num: %d' % (msa2, m2.seqlen, m2.seqNum)

	print 'loading map ..'
	pairdict = {}
	count = 0
	with open(pairfile) as fp:
		for line in fp:
			strArray = line.split(' ')
			# format: primary_name suffix_name
			pairdict[strArray[0]] = strArray[1].strip()
			count+=1
	print repr(pairdict)
	print '%d map loaded.' % count

	count=0
	fout = open(outfile, 'w')
	for mi in m1.msaArray:
		hi = mi[0] # >sp|Q14126|DSG2_HUMAN Desmoglein-2 OS=Homo sapiens GN=DSG2 PE=1 SV=2
		hArrayi = hi.split('|') # Q14126

		prefix=hArrayi[1]
		suffix = pairdict[prefix]
		si = mi[1]

		for mj in m2.msaArray:
			hj = mj[0]
			hArrayj = hj.split('|') # Q14126

			sj = mj[1]
			if suffix == hArrayj[1]:
				header = '>%s|%s' % (prefix, suffix)
				seq = '%s%s' % (si,sj)
				fout.write('%s\n%s\n' % (header, seq))
				count+=0

	print '%d pairs have been combined.' % count
	fout.close()

if __name__ == '__main__':
	main()