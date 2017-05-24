import sys

def sliceHTML(line, token):
	pos1 = line.find(token)
	if pos1!= -1:
		pos2 = line.find('"',pos1+len(token))
		if pos2 == -1:
			print 'error::umbiguous line:\n%s for token [%s]' % (line, token)
			exit(-1)
		return line[pos1+len(token):pos2]
	else:
		return False

def extractpfamid(arglist):
	if len(arglist) < 2:
		print 'Usage: python proc_extractPfamID.py 1aos.rcsb'
		print 'output: 1oas.rcsb.pfam'
		return

	pdbname = arglist[2][:-5] 
	ctoken = 'tr class="chain'
	ptoken = 'pfamID='
	start = 'Protein Family Annotation'
	end = 'Gene Product Annotation'
	chain = []
	pfam = []
	flag = False
	with open(arglist[2]) as fp:
		for line in fp:
			if start in line:
				flag = True
			if flag == False:
				continue

			if end in line:
				break

			c = sliceHTML(line, ctoken)
			if c!=False:
				chain.append(c)
				continue
			p = sliceHTML(line, ptoken)
			if p!=False:
				pfam.append(p)
				continue
	if len(pfam) != len(chain):
		print 'error::unmatched pair'
		print repr(pfam)
		print repr(chain)
		return

	for i in xrange(0, len(pfam)):
		print '%s %s %s' % (pdbname, chain[i], pfam[i])


def main():
	extractpfamid([' ']+sys.argv)

if __name__ == '__main__':
	main()