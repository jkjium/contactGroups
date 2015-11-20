import sys
from sparseIsing import sparseIsing
def main():
	#if len(sys.argv) < 2:
	#	print 'Usage python proc_spdb.py XXXX_A.domain'
	#	return
	si = sparseIsing(2.7, 4)
	print 'lambda: %f' % si.la
	print 'a: %f' % si.a
	print 't: %f' % 2
	print 'I_le_la: %f' % si.I_le_la(2)
	print 'I_g_la: %f' % si.I_g_la(2)
	print 'p_lambda_t: %f' % si.P_lambda_t(2)
	#si.W_()
	print si.beta
	#print si.w
	print 'S_r_t: %f' % si.S_r_t(-3, 2)
	print si.idict

if __name__=="__main__":
	main()