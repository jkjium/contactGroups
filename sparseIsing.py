import sys
import numpy as np

class sparseIsing(object):

    def __init__(self, la, N):
    	self.N = N # 20A + 20AA + 20AAA
    	self.la = la # lambda 
    	self.a = 3.7 # p3 Following Fan and Li (2001) set a=3.7

    	self.beta =  np.random.uniform(-3.0, 3.0, self.N * (self.N - 1) / 2)



    def I_le_la(self, t):
    	if t <= self.la:
    		return 1
    	else:
    		return 0

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

    def P_lambda_(self,t):
    	return self.la * (self.I_le_la(t) + self.tp(self.a * self.la - t) / ((self.a - 1) * self.la) * self.I_g_la(t))


