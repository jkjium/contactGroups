import sys
from utils_pfammsa import utils_pfammsa

# updated from proc_extractMSAbyname.py
def extractbyname(arglist):
	if len(arglist)!=5:
		print 'Usage: python proc_pfammsa.py extractp90 PF00008_full.txt p90_PF00008.title PF00008_p90.txt'
		return

	msafile = sys.argv[2]
	titlefile = sys.argv[3]
	outfile = sys.argv[4]

	count=0
	up = utils_pfammsa(msafile)
	fout = open(outfile ,'w')
	with open(titlefile) as fp:
		for line in fp: # >A0A074XN79_AURPU/1192-1206 A0A074XN79.1 PF02809.19;UIM;
			strarr = line.split(' ')
			entry = up.fetchbyname(strarr[0][1:])
			if entry!=False:
				fout.write(entry)
				count+=1
	fout.close()
	print 'save %d seqs in %s.' % (count, outfile)


def main():
	dispatch = {
		'extractp90': extractbyname
	}

	cmd = sys.argv[1]
	if cmd not in dispatch:
		print 'Invalid cmd %s' % cmd
		return

	dispatch[sys.argv[1]](sys.argv)	

if __name__ == '__main__':
	main()
