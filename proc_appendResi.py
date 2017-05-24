import sys
from collections import defaultdict
"""
For fdr project
append active site resi to categorized pdb pfam map file

0-CSA_2_0_121113.txt:
13pk,1,Arg,B,39,S,LIT,13pk

1-cas_pdb_chain.all.pfam:
132l A PF00062
135l A PF00062

"""

def main():
	if len(sys.argv) < 4:
		print 'Usage: python proc_appendResi.py 0-CSA_2_0_121113.txt 1.5-cas_pdb_chain.all.pfam 2-pdb_chain_pfam_resi.txt'
		return

	outfile = sys.argv[3]


	fulldict = defaultdict(lambda:[])

	# read full list
	# 13pk,0,Arg,A,39,S,LIT,13pk
	with open(sys.argv[1]) as fp:
		for line in fp:
			strarr = line.strip().split(',')
			key = '%s %s' % (strarr[0], strarr[3]) # pdb + ' ' + chain
			value = '%s%s' % (strarr[2], strarr[4]) # resn + resi
			fulldict[key].append(value)

	print '%d unique pdb+chains collected.' % len(fulldict)
	'''
	for k in fulldict:
		print '%s: %s' % (k, repr(fulldict[k]))
	'''

	# read 1-cas_pdb ...
	# 12as A PF03590
	count = 0
	fout = open(sys.argv[3], 'w')
	with open(sys.argv[2]) as fp:
		for line in fp:
			strarr = line.strip().split(' ')
			key = '%s %s' % (strarr[0], strarr[1])
			if len(fulldict[key])==0:
				print 'error::%s resi not found.' % key
				continue
			outstr = '%s %s\n' % (line.strip(), ','.join(fulldict[key]))
			fout.write(outstr)
			count+=1
	fout.close()

	print '%d records append.' % count
	print 'save to %s.' % outfile

if __name__ == '__main__':
	main()