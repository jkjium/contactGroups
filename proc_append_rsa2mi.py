import sys
import commp as cp

def main():
	if len(sys.argv) < 4:
		cp._err('Usage: proc_append_rsa2mi.py rsafile mifile outfile')

	rsafile = sys.argv[1]
	mifile = sys.argv[2]
	outfile = sys.argv[3]

	rsadict = {}
	with open(rsafile) as fp:
		for line in fp:
			line = line.strip()
			if len(line) == 0:
				continue
			# TYR A 117 56.22 26.4
			sarr = line.split(' ')
			rsadict[sarr[2]] = sarr[4]
	cp._info('load %d entries' % len(rsadict))

	outlist = []
	with open(mifile) as fp:
		for line in fp:
			line = line.strip()
			if len(line) == 0:
				continue
			# 171 172 119 120 1.62
			sarr = line.split(' ')
			outlist.append('%s %s %s' % (line, rsadict[sarr[0]], rsadict[sarr[1]]))
	cp._info('map %d entries' % len(outlist))

	with open(outfile, 'w') as fout:
		fout.write('%s' % ('\n'.join(outlist)))
	cp._info('save to %s' % outfile)

if __name__ == '__main__':
	main()