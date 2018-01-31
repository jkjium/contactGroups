import sys
import commp as cp

def main():
	if len(sys.argv)!=5:
		cp._err('Usage: python proc_batch_gen.py dbname stubfile sm_name flag{preset|alter}')
	# stubfile contains all the .fa filenames
	dbname = sys.argv[1]
	stubfile = sys.argv[2]
	sm = sys.argv[3]
	flag = sys.argv[4]

	if flag not in ['preset', 'alter']:
		cp._err('Usage: python proc_batch_gen.py dbname stubfile sm_name flag{preset|alter}')

	# $ blastp -query t.fa -db $ASTRALS40 -outfmt "10 stitle evalue" -evalue 0.0001 -matrix BLOSUM62 -gapopen 9 -gapextend 2 -out t.out
	fout = open('batch_blast.sh','w')
	with open(stubfile) as fp:
		for line in fp:
			line = line.strip()
			if len(line)==0:
				continue
			#for g in gap:
			for g in cp.gapb80:
				#fout.write('blastp -query %s -db $%s -outfmt "10 stitle evalue" -evalue 0.01 -matrix BLOSUM62 -gapopen %d -gapextend %d -out %s.%s.%d-%d.out\n' % (line, dbname, g[0], g[1], line, sm, g[0], g[1]))
				if flag == 	'preset':
					fout.write('blastp -query %s -db $%s -outfmt "10 stitle evalue" -evalue 10 -matrix %s -gapopen %d -gapextend %d -out %s.%s.%d-%d.out\n' % (line, dbname, sm, g[0], g[1], line, sm, g[0], g[1]))
				else flag == 'alter':
					fout.write('blastp -query %s -db $%s -outfmt "10 stitle evalue" -evalue 10 -matrix BLOSUM62 -gapopen %d -gapextend %d -out %s.%s.%d-%d.out\n' % (line, dbname, g[0], g[1], line, sm, g[0], g[1]))
	fout.close()
	cp._info('batch_blast.sh generated for %s' % sm)

if __name__ == '__main__':
	main()
