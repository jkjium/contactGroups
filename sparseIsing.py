import sys
import numpy as np

class sparseIsing(object):

    def __init__(self, la, N):
    	self.N = N # 20A + 20AA + 20AAA
    	self.la = la # lambda 
    	self.a = 3.7 # p3 Following Fan and Li (2001) set a=3.7

    	self.beta =  np.random.uniform(-5.0, 5.0, self.N * (self.N - 1) / 2)
        self.w = []

        self.idict= {}
        self.i2jk()

    # convert linear index to pairwised indices
    def i2jk(self):
        count = 0
        for i in xrange(1, self.N+1):
            for j in xrange(i+1, self.N+1):
                self.idict[count] = (i, j)
                count+=1


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


    # p5 S(r,t) = sgn(r)(|r|-t)+
    def S_r_t(self, r, t):
        return np.sign(r) * (self.tp(abs(r) - t))

    # sequentially get new beta vector by S_r_t
    def S_beta(self):
        for i in xrange(0, len(self.beta)):
            r = self.beta[i] + 2 * self.Z_jk()
            t = 2 * self.P_lambda_t(abs(self.beta[i]))
            self.beta[i] = self.S_r_t(r, t)

        return 



