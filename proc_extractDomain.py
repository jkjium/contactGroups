# pre-process all atom pdb
# extract CA in Chain A and save to file
import sys
from protein import protein
from atom import atom
def main():
	if len(sys.argv) < 2:
		print "Usage python proc_extractDomain.py domain_desc_file"
		return

	fin = open(sys.argv[1], 'r')
	for line in fin.readlines():
		line = line.strip()
		strArr = line.split(',')

		tip_filename = strArr[0]
		start = int(strArr[1])
		end = int(strArr[2])

		fo = open(tip_filename+'.domain', 'w')
		p = protein(tip_filename, 'alpha',center='TIP') 

		for a in p.atoms:
			if a.resSeq >= start && a.resSeq <= end:
				fo.write(a.writeAtom())

		fo.close()

	fin.close()
if __name__=="__main__":
	main()

