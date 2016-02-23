'''
Extract naccess information
work with cgroup class
rsa format:

REM  Relative accessibilites read from external file "standard.data"
REM  File of summed (Sum) and % (per.) accessibilities for 
REM RES _ NUM      All-atoms   Total-Side   Main-Chain    Non-polar    All polar
REM                ABS   REL    ABS   REL    ABS   REL    ABS   REL    ABS   REL
RES ILE B 397    96.91  55.3  59.11  42.8  37.80 101.7  59.76  42.9  37.15 103.3
RES ASP B 398     5.00   3.6   0.19   0.2   4.82  12.8   0.27   0.5   4.73   5.2
...
...
RES ASP B 653   111.79  79.6  98.52  95.9  13.28  35.2  51.17 103.9  60.63  66.5
RES VAL B 654    80.44  53.1  12.88  11.3  67.56 181.8  41.41  35.9  39.02 108.5
END  Absolute sums over single chains surface 
CHAIN  1 B    13134.9      11364.2       1770.7       7723.4       5411.4
END  Absolute sums over all chains 
TOTAL         13134.9      11364.2       1770.7       7723.4       5411.4

'''

import numpy as np
import time
import sys
from AAmap import AAmap
from collections import defaultdict

__all__=['naccess', 'rsa']

'''
class that describe a single line in rsa file
'''
class rsa(object):
	def __init__(self, line):
		strArr = line.split()
		self.resn = strArr[1]
		self.chain = strArr[2]
		self.resi = strArr[3]

		# All-atoms
		self.AA_ABS = float(strArr[4])
		self.AA_REL = float(strArr[5])

		# Total-Side
		self.TS_ABS = float(strArr[6])
		self.TS_REL = float(strArr[7])

		# Main-Chain
		self.MC_ABS = float(strArr[8])
		self.MC_REL = float(strArr[9])

		# Non-polar
		self.NP_ABS = float(strArr[10])
		self.NP_REL = float(strArr[11])

		# All-polar
		self.AP_ABS = float(strArr[12])
		self.AP_REL = float(strArr[13])

	# dump content string
	def outString(self):
		return '%s - %s - %s - %.2f - %.2f - %.2f - %.2f - %.2f - %.2f - %.2f - %.2f - %.2f - %.2f' % \
				(self.resn, self.chain, self.resi, 
					self.AA_ABS, self.AA_REL, 
					self.TS_ABS, self.TS_REL, 
					self.MC_ABS, self.MC_REL, 
					self.NP_ABS, self.NP_REL, 
					self.AP_ABS, self.AP_REL)



'''
class that store a whole ras file information
'''
class naccess(object):

	# nafile name must in the format of "XXXX.rsa"
	def __init__(self, nafile):
		self.pdb = nafile[0:4]
		self.rsaDict = {}
		self.resiDict = defaultdict(lambda: '')
		self.alphabet = ['B', 'E']
		aamap = AAmap()

		lines = [line.strip() for line in open(nafile)]
		for naline in lines:
			head = naline[0:3]
			if head == 'RES':
				r = rsa(naline)
				key = '%s%s%s' % (aamap.getAAmap(r.resn), r.chain, r.resi)
				self.rsaDict[key] = r

				varkey = '%s%s' % (aamap.getAAmap(r.resn), self.accessible(key))
				varvalue = '%s%s%s ' % (self.resiDict[varkey], r.chain, r.resi)
				self.resiDict[varkey] = varvalue
			elif head == 'TOTAL':
				key = 'TOTAL'
				self.rasDict[key] = naline.split()


	# dump naccess content
	def dump(self):
		for key in self.rsaDict:
			print '[%s]: %s' % (key, self.rsaDict[key].outString())


	# alphabet assigning
	def accessible(self, key):
		if self.rsaDict[key].AA_REL > 20.0:
			return 'E'
		else:
			return 'B'

	# output residue number according to query variable type
	def dumpResiMap(self):
		for key in self.resiDict:
			print '[%s]: %s' % (key, self.resiDict[key])


	# output resi list for varname
	# example: input 'DB'
	# return: [DB]: B398 B504 B521 B579
	def getResiList(self, varname):
		return self.resiDict[varname]