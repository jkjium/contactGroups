import sys
import commp as cp

def searchbytitle(msafile, title):
	for head, seq in cp.fasta_iter(msafile):
		if head == title:
			return seq
	return False
"""
get seq from seed and p90 MSA by referring common title
output the map from p90 position to seed position
"""
def main():
	if len(sys.argv)!=2:
		print 'Usage: python proc_p902seedmap.py PF00008'
		print 'input: PF00008.comm, PF00008_p90.txt, PF00008.seed'
		print 'output: map: p90 pos -> seed pos'
		return

	prefix = sys.argv[1]
	commfile = prefix + '.comm'
	p90file = prefix + '_p90.txt'
	seedfile = prefix + '.seed'
	outfile = prefix + '_p902seed.posmap'

	with open(commfile) as fp:
		titles = fp.readlines()

	posmap = {}
	for line in titles:
		title = line[1:].strip()
		seedseq = searchbytitle(seedfile, title)
		p90seq = searchbytitle(p90file, title)
		if seedseq==False or p90seq==False:
			continue
		posmap = cp.posmap(p90seq.upper(), seedseq.upper())
		if posmap!=False:
			break

	if len(posmap)==0:
		print 'error: %s no title matched.' % prefix
		return

	fout = open(outfile, 'w')
	for k in posmap:
		fout.write('%d %d\n' % (k, posmap[k]))
	fout.close()
	print 'save file: %s' % outfile

if __name__ == '__main__':
	main()