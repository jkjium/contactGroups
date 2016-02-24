'''
get msa position id 
'''
import sys
from msa import msa
from protein import protein

def main():
	if len(sys.argv)<4:
		print 'Usage: utils_naccess.py msafile pdbfile chain+resi'
		return
	msafile = sys.argv[1]
	pdbfile = sys.argv[2]
	resi = sys.argv[3]

	m = msa(msafile)
	p = protein(pdbfile)
	#print p.seq
	#print p.resDict
	resIdx = p.resDict[resi]
	posMap = m.getPosMap(p)
	#print m.msaArray[0][1]
	for i in posMap:
		if i == resIdx[0]:
			print '[Res: %s] : (seqi: %d (%s) - msai: %d (%s))' % (resi, resIdx[0], resIdx[1], posMap[i], m.msaArray[0][1][posMap[i]])
			break

if __name__ == '__main__':
	main()
