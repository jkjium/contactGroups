# extract cath domain from tip file by referening domain_desc_file (pdbname, startRES, endRES)
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
		
		print tip_filename+'.domain'
		fo = open(tip_filename+'.domain', 'w')
		p = protein('a'+tip_filename+'.tip', 'alpha',center='TIP') 

		for a in p.atoms:
			if a.resSeq >= start and a.resSeq <= end:
				fo.write(a.writeAtom())

		fo.close()

	fin.close()
if __name__=="__main__":
	main()

