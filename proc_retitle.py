import sys
import commp as cp

def main():
	if len(sys.argv)!=4:
		cp._er('Usage: proc_retitle.py orig.fa domain.list out.fa')
	
	origfafile = sys.argv[1]
	domainfile = sys.argv[2]
	outfile = sys.argv[3]

	domaindict = {}
	# load domain dictionary
	with open(domainfile) as fp:
		for line in fp:
			line = line.strip()
			if len(line) ==0:
				continue
			sarr = line.split()
			domaindict[sarr[0]] = '%s-%s-%s-%s' % (sarr[1], sarr[2], sarr[3], sarr[4])
	#print repr(domaindict)

	# iterate the original fa sequences and replace the title
	outlist = []
	'''
	$ grep ">" cath-dataset-nonredundant-S20.fa|awk '{print substr($0, 0, 12)}'|sort|uniq -c
  	14433 >cath|4_2_0|	
	header: >cath|4_2_0|12asA00/4-330
	'''
	fout = open(outfile, 'w')
	for header, seq in cp.fasta_iter(origfafile):
		#print header[11:18], domaindict[header[11:18]]
		fout.write('>%s %s\n%s\n' % (header[11:18], domaindict[header[11:18]], seq))
	cp._info('save to %s.' %  outfile)
	fout.close()
		

if __name__=="__main__":
	main()
