import sys
import numpy as np

class sparseIsing(object):

    def __init__(self, la, dataFile):
    	self.dataFile = dataFile
        self.data = {} # boolean value for P variables with key = {0, 1, .... N}
    	self.N, self.P = self.loadData(dataFile) 
    	self.la = la # lambda 
    	#self.a = 3.7 # p3 Following Fan and Li (2001) set a=3.7
    	self.a = 3.7 # p3 Following Fan and Li (2001) set a=3.7

    	self.beta =  np.random.uniform(-5.0, 5.0, self.P * (self.P - 1) / 2)
        self.w = np.empty_like(self.beta)

        # index map between linear index and row, column index
        self.i2jk = {}
        self.jk2i = {}
        self.init_indexMap()


    def loadData(self, filename):
    	fp = open(filename, 'r')
    	lines = fp.readlines()
    	fp.close()
    	i = 0
    	for line in lines:
    		strArr = line.strip().split(' ')
    		self.data[i] = [int(k) for k in strArr]
    		i+=1
    	p = len(strArr)
    	return i, p

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
    	print 'r:%f, t:%f, S_r_t:%f' % (r,t,self.S_r_t(r,t))
    	return self.S_r_t(r, t)


    # p6 algorithm 1. sequentially update beta_ij
    def CMA(self):
    	count=0
    	# iterate until converge
    	while True:
    		beta_prev = np.empty_like(self.beta)
    		beta_prev[:] = self.beta
    		for i in xrange(0, len(self.beta)):
    			(j, k) = self.i2jk[i]
    			self.beta[i] = self.S_beta_j_k(j, k)
    		d = np.linalg.norm(self.beta - beta_prev)
    		print 'CMA::iter[%d]::d: %f' % (count, d)
    		count+=1
    		if d < 1e-4:
   				#print self.beta
  				return

  	# p8 main algorithm LLA-CMA		
    def LLA_CMA(self):
    	count = 0
    	# initialize w_jk
    	for i in xrange(0, len(self.w)):
    		self.w[i] = self.P_lambda_t(abs(self.beta[i]))

    	# iterate step
    	while True:
    		w_prev = np.empty_like(self.w)
    		w_prev[:] = self.w
    		self.CMA() # update all beta
    		for i in xrange(0, len(self.w)):
    			self.w[i] = self.P_lambda_t(abs(self.beta[i]))
    		d = np.linalg.norm(self.w - w_prev)
    		print 'LLA::iter[%d]::d: %f' % (count, d)
    		count+=1
    		if d < 1e-4:
    			#print self.w
    			return

    # output result
    def formatResult(self):
    	fp = open(self.dataFile+'.result', 'w')
    	for i in xrange(0, len(self.w)):
    		(j, k) = self.i2jk[i]
    		fp.write('beta(%d, %d): %f w(%d, %d): %f\n' % (j, k, self.beta[i], j, k, self.w[i])) 
    		print 'beta(%d, %d): %f w(%d, %d): %f\n' % (j, k, self.beta[i], j, k, self.w[i]) 
