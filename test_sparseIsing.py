import sys
from sparseIsing import sparseIsing
def main():
	#if len(sys.argv) < 2:
	#	print 'Usage python proc_spdb.py XXXX_A.domain'
	#	return
	si = sparseIsing(2.7, 4, 'test_sparseIsing.data')
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
	print si.i2jk
	print si.jk2i
	print si.data
	print si.N
	print si.theta_j_n(0,1)
	print si.z_j_k(0,1)
	print si.S_beta_j_k(0,1)

if __name__=="__main__":
	main()