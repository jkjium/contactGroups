import sys
import commp  as cp

def main():
	if len(sys.argv) < 4:
		cp._err('Usage: python proc_cathfafilter.py clffile fafile outfile')
	clffile = sys.argv[1]
	fafile = sys.argv[2]
	outfile = sys.argv[3]

	# load description file
	cp._info('loading clf file ...')
	with open(clffile) as fp:
		# 1oaiA00     1    10     8    10     1     1     1     1     1    59 1.000
		'''
		Column 1:  CATH domain name (seven characters)
		Column 2:  Class number
		Column 3:  Architecture number
		Column 4:  Topology number
		Column 5:  Homologous superfamily number
		Column 6:  S35 sequence cluster number
		Column 7:  S60 sequence cluster number
		Column 8:  S95 sequence cluster number
		Column 9:  S100 sequence cluster number
		Column 10: S100 sequence count number
		Column 11: Domain length
		Column 12: Structure resolution (Angstroms)
		           (999.000 for NMR structures and 1000.000 for obsolete PDB entries)
		'''
		cathdict = {}
		for line in fp:
			line = line.strip()
			if len(line)==0 or line[0] == '#':
				continue
			sarr = line.split()
			name = sarr[0]
			cathid = '-'.join(sarr[1:5])
			cathdict[name] = cathid
	cp._info('%d records loaded.' % len(cathdict))

	# iterate fa file
	fout = open(outfile, 'w')
	cp._info('filtering fa ...')
	for head, seq in cp.fasta_iter(fafile):
		#cath|current|1o6aA00/68-154
		sarr = head.split('|')
		name = sarr[2][0:7]
		if name not in cathdict:
			cp._err('no desc found for %s' % name)
		fout.write('>%s %s\n%s\n' % (name, cathdict[name], seq))
	fout.close()
	cp._info('save to %s' % outfile)

if __name__ == '__main__':
	main()