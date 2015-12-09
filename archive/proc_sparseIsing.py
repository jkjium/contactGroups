import sys
from sparseIsing import sparseIsing
def main():
	if len(sys.argv) < 2:
		print 'Usage python proc_sparseIsing.py matrixFile'
		return
	inFile = sys.argv[1]
	si = sparseIsing(2.7, inFile)
	print 'inFile: %s' % inFile
	print 'lambda: %f' % si.la
	print 'a: %f' % si.a
	print 'N: %d' % si.N
	print 'P: %d' % si.P
	si.LLA_CMA()
	si.formatResult()

if __name__=="__main__":
	main()