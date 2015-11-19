import sys
from sparseIsing import sparseIsing
def main():
	#if len(sys.argv) < 2:
	#	print 'Usage python proc_spdb.py XXXX_A.domain'
	#	return
	si = sparseIsing(2.7)
	print 'lambda: %f' % si.la
	print 'a: %f' % si.a
	print 't: %f' % 2
	print 'I_le_la: %f' % si.I_le_la(2)
	print 'I_g_la: %f' % si.I_g_la(2)
	print 'p_lambda_: %f' % si.P_lambda_(2)

if __name__=="__main__":
	main()