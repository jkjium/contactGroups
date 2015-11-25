import sys
import numpy as np

class sparseIsing(object):

    def __init__(self, la, P, dataFile):
    	self.P = P # 20A + 20AA + 20AAA
    	self.la = la # lambda 
    	self.a = 3.7 # p3 Following Fan and Li (2001) set a=3.7

    	self.beta =  np.random.uniform(-5.0, 5.0, self.P * (self.P - 1) / 2)
        self.w = []

        # index map between linear index and row, column index
        self.i2jk = {}
        self.jk2i = {}
        self.init_indexMap()

        self.data = {} # boolean value for P variables
    	self.N = self.loadData(dataFile) 

    def loadData(self, filename):
    	fp = open(filename, 'r')
    	lines = fp.readlines()
    	fp.close()
    	i = 0
    	for line in lines:
    		strArr = line.strip().split(' ')
    		self.data[i] = [int(k) for k in strArr]
    		i+=1
    	return i

    # convert linear index to pairwised indices
    def init_indexMap(self):
        i = 0
        for j in xrange(0, self.P):
            for k in xrange(j+1, self.P):
                self.i2jk[i] = (j, k)
                self.jk2i['%d,%d' % (j,k)] = i
                self.jk2i['%d,%d' % (k,j)] = i
                i+=1

    # I(t<=L)
    def I_le_la(self, t):
    	if t <= self.la:
    		return 1
    	else:
    		return 0

    # I(t>L)
    def I_g_la(self, t):
    	if t > self.la:
    		return 1
    	else:
    		return 0

    # truncate positive
    def tp(self, v):
    	if v > 0:
    		return v
    	else:
    		return 0

    # calculate P_lambda_t for single value
    def P_lambda_t(self,t):
    	return self.la * (self.I_le_la(t) + self.tp(self.a * self.la - t) / ((self.a - 1) * self.la) * self.I_g_la(t))

    # get Wij vector
    #def W_(self):
    #    for i in xrange(0, len(self.beta)):
    #        self.w.append(self.P_lambda_t(self.beta[i]))

    # p2 theta_j_n; jth element and nth observation
    def theta_j_n(self, j, n):
    	v = 0.0
    	d = self.data[n]
    	for k in xrange(0, self.P):
    		if j != k:
				x_jn = d[j]
				x_kn = d[k]
				beta_jk = self.beta[self.jk2i['%d,%d' % (j,k)]]
				v+=beta_jk*x_jn*x_kn
		return np.exp(v)/(np.exp(v)+1)


	# p5 (2.4) z_j_k
    def z_j_k(self, j, k):
    	v = 0.0
    	for n in xrange(0, self.N):
    		d = self.data[n]
    		x_jn = d[j]
    		x_kn = d[k]
    		theta_jn = self.theta_j_n(j, n)
    		theta_kn = self.theta_j_n(k, n)
    		v+=x_kn*x_jn*(2-theta_kn-theta_jn)
    	return v/self.N
	

    # p5 S(r,t) = sgn(r)(|r|-t)+
    def S_r_t(self, r, t):
        return np.sign(r) * (self.tp(abs(r) - t))


    # p6 algorithm 1 step (2)
    # sequentially get new beta vector by S_r_t
    def S_beta_j_k(self, j, k):
    	i = self.jk2i['%d,%d' % (j,k)]
    	r = self.beta[i] + 2 * self.z_j_k(j, k)
    	t = 2 * self.P_lambda_t(abs(self.beta[i]))
    	print 'r:%f, t:%f' % (r,t)
    	return self.S_r_t(r, t)

    #    for i in xrange(0, len(self.beta)):
    #        r = self.beta[i] + 2 * self.Z_jk()
    #        t = 2 * self.P_lambda_t(abs(self.beta[i]))
    #        self.beta[i] = self.S_r_t(r, t)

    #    return 



