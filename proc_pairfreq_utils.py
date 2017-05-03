import sys
import numpy as np
import operator
import math
import common.commp as cm
from collections import defaultdict


def normalize(k,v,p):
	'''
	read .pairfreq list and sum up all the items
	called in pairfreqsum()
	'''
	pk1 = ''.join(sorted([k[0], k[1]]))
	pk2 = ''.join(sorted([k[2], k[3]]))
	return float(v)/(p[pk1]+p[pk2])

def pcount(psm):
	'''
	calculate substitution count for residue contact
	called in pairfreqsum
	'''
	p1 = defaultdict(lambda:0) # single changed case
	p2 = defaultdict(lambda:0) # double changed case

	for k in psm:
		#if k[0]=='.' and k[1]=='.':
		#	continue
		#elif (k[0]!=k[2]) and (k[1]!=k[3]):
		if (k[0]!=k[2]) and (k[1]!=k[3]):
			# for both changed case
			# DACG - p2k1: AD - p2k2: CG
			p2k1 = ''.join([k[0], k[1]])
			#p2k1 = ''.join(sorted([k[0], k[1]]))
			p2[p2k1]+=psm[k]

			#p2k2 = ''.join([k[2], k[3]])
			p2k2 = ''.join(sorted([k[2], k[3]]))
			p2[p2k2]+=psm[k]

		elif (k[0]==k[2] and (k[1]==k[3])):
			# does not need to normalize
			continue
		else:
			# single changed case
			p1k1 = ''.join([k[0], k[1]])
			#p1k1 = ''.join(sorted([k[0], k[1]]))
			p1[p1k1]+=psm[k]

			p1k2 = ''.join([k[2], k[3]])
			#p1k2 = ''.join(sorted([k[2], k[3]]))
			p1[p1k2]+=psm[k]

	return p1,p2


def pairfreqsum():
	"""
	split .pairfreq into .pairfreq.2 , .1, .0
	"""
	if len(sys.argv)<3:
		print 'Usage: python proc_pairfreq_utils.py pairfreqsum pairfreq.list'
		return

	listfile = sys.argv[2]
	outfile = listfile+'.allpairfreq'

	psm = defaultdict(lambda:0)
	with open(listfile) as fp:
		for pfile in fp:
			pfile = pfile.strip()
			if len(pfile) < 2:
				continue
			with open(pfile) as fp1:
				for line in fp1:
					strarr = line.strip().split(' ')
					if len(strarr) < 2 or len(set(strarr[0]).intersection(cm.abaa))!=0:
						continue
					lowerlist = [c for c in strarr[0] if c.islower()]
					if len(lowerlist)!=0:
						continue
					psm[strarr[0]]+=int(strarr[1])
			print '%s processed.' % pfile

	print 'calculating normalization terms ...'
	p1,p2 = pcount(psm)

	fo0 = open(outfile+'.0', 'w')
	fo1 = open(outfile+'.1', 'w')
	fo2 = open(outfile+'.2', 'w')

	print 'save to %s.0' % outfile
	print 'save to %s.1' % outfile
	print 'save to %s.2' % outfile

	for k,v in sorted(psm.items(), key=operator.itemgetter(0), reverse=True): # sort by key
		if k[0]=='.' and k[1]=='.':
			continue
		if (k[0]!=k[2]) and (k[1]!=k[3]):
			outstr = '%s %d %f\n' % (k, v, normalize(k,v,p2))
			fo2.write(outstr)
		elif (k[0]==k[2] and (k[1]==k[3])):
			outstr = '%s %d\n' % (k, v)
			fo0.write(outstr)
		else:
			outstr = '%s %d %f\n' % (k, v, normalize(k,v,p1))
			fo1.write(outstr)

	fo0.close()
	fo1.close()
	fo2.close()


def sentropy(d):
	"""
	calculate information entropy for the given dictionary
	called in ssentropy()
	"""
	#sort by key
	values = [v for k,v in sorted(d.items(), key=operator.itemgetter(0))]
	total = sum(values)
	if total == 0.0:
		return 0.0
	plist = np.array(values)/float(total)

	return -sum([p*np.log2(p) for p in plist if p!=0.0])

def ssentropy():
	"""
	calculate the neighbor 20x(20-1)/2=190 aa composition and the entropy for one type of AA
	"""
	if len(sys.argv)<3:
		print 'Usage: python proc_pairfreq_utils.py ssentropy pairfreq.1.top'
		print 'output: pairfreq.1.top.ssentropy'
		return

	infile = sys.argv[2]
	outfile = infile + '.ssentropy'

	smdict = {}
	print 'loading singled subsitutions ...'
	with open(infile) as fp:
		for line in fp:
			# GPNP 12882049 7.1100
			strarr = line.strip().split(' ')
			if len(strarr)!=3:
				print 'error::incomplete line [%s]' % line
				continue
			if len(set(strarr[0]).intersection(cm.abaa))!=0:
				#print 'info::gap [%s] found. skip.' % strarr[0]
				continue
			# get the shared character as the mainkey
			# sorted substitution as the subkey
			strar = strarr[0]
			if strar[0] == strar[2]:
 				mainkey = strar[0]
 				subkey = ''.join(sorted([strar[1], strar[3]]))
 			else:
 				mainkey = strar[1]
 				subkey = ''.join(sorted([strar[0], strar[2]]))

 			# add up the normalized substitution frequency
 			if mainkey not in smdict:
 				d = defaultdict(lambda:0) # init 190 possible neighbors
 				'''
 				for i in xrange(0, len(cm.aas01)):
 					for j in xrange(i+1, len(cm.aas01)):
 						a = '%s%s' % (cm.aas01[i], cm.aas01[j])
 						d[a]=0.0
 				'''
 				d[subkey] = float(strarr[2])
 				smdict[mainkey] = d
 			else:
 				t = smdict[mainkey]
 				t[subkey]+=float(strarr[2])
 				#smdict[mainkey][subkey]+=float(strarr[2])

 	# iterate again to calculate the total and the entropy
 	count = 0
 	fout = open(outfile, 'w')
 	for mk in smdict:
 		sd = smdict[mk] # sd is the dictionary for 190 neighboring substitution
 		h = sentropy(sd)
 		total = sum([sd[k] for k in sd])
 		strn190 = ' '.join(['%s %f' % (k, v) for k,v in sorted(sd.items(), key=operator.itemgetter(1), reverse=True)])
 		outstr = '%s %f %f %s\n' % (mk, total, h, strn190)
 		fout.write(outstr)
 		count+=1
 		print '%d/%d [%s]' % (count, len(smdict), mk)
 	print 'save to [%s].' % outfile
 	fout.close()



def composition():
	"""
	bundle all the substitutions with the same amino acid compositions
	"""
	if len(sys.argv)< 3:
		print 'Usage: python proc_pairfreq_utils.py composition pairfreq.2.top'
		print 'output: pairfreq.2(1).top.composition'
		return

	infile =sys.argv[2]
	outfile = infile+'.composition'

	smdict = {}
	with open(infile) as fp:
		for line in fp:
			strarr = line.strip().split(' ')
			if len(strarr) != 3:
				print 'error::incomplete line [%s]' % line
				continue
			if len(set(strarr[0]).intersection(cm.abaa))!=0:
				continue
			s = sorted(strarr[0])
			mainkey = ''.join(s)
			if mainkey not in smdict:
				# initialize with all the combinations (last three position)
				smdict[mainkey] = defaultdict(lambda:[0,0])
				'''
				smdict[mainkey][''.join([s[0],s[1],s[2],s[3]])]=[0,0]
				smdict[mainkey][''.join([s[0],s[1],s[3],s[2]])]=[0,0]
				smdict[mainkey][''.join([s[0],s[2],s[1],s[3]])]=[0,0]
				smdict[mainkey][''.join([s[0],s[2],s[3],s[1]])]=[0,0]
				smdict[mainkey][''.join([s[0],s[3],s[1],s[2]])]=[0,0]
				smdict[mainkey][''.join([s[0],s[3],s[2],s[1]])]=[0,0]
				'''

				# with the minimum one be at the first position
				# here should not have non-existed key
				smdict[mainkey][strarr[0]] = [int(strarr[1]), float(strarr[2])]
			else:
				smdict[mainkey][strarr[0]][0]+=int(strarr[1])
				smdict[mainkey][strarr[0]][1]+=float(strarr[2])

	fout = open(outfile, 'w')
	count = 0
	for mk in smdict:
		ps = smdict[mk]
		countsum = sum([ps[k][0] for k in ps])
		normsum = sum([ps[k][1] for k in ps])

		outstr = '%f %d %s\n' % (normsum, countsum, ' '.join(['%s %d %f' % (k,ps[k][0],ps[k][1]) for k in ps]))
		fout.write(outstr)
		count+=1
		print '%d/%d' % (count, len(smdict))

	print 'save to %s.' % outfile
	fout.close()

def dist210():
	"""
	calculate pair substitution distribution across all 20x(20-1)/2 possible pairs
	"""
	if len(sys.argv)<3:
		print 'Usage: python proc_pairfreq_utils.py dist210 pairfreq.list.allpairfreq'
		print 'Output: pairfreq.list.allpairfreq.dist210'
		return

	infile = sys.argv[2]
	outfile = infile+'.dist210'

	smdict = {}
	print 'dispatching ..'
	with open(infile) as pf:
		# WYYW 155254 0.000418
		for line in pf:
			strarr = line.strip().split(' ')
			if len(strarr) != 3:
				continue

			subdict = {}
			for i in xrange(0, len(cm.aas01)):
				for j in xrange(0, len(cm.aas01)):			
					subdict['%s%s' % (cm.aas01[i], cm.aas01[j])] = 0.0

			mainkey = strarr[0][0:2] # should be sorted
			subkey = strarr[0][2:]
			rawcount = float(strarr[1]) # use raw count instead of the nfreq

			if mainkey not in smdict:
				subdict[subkey] = rawcount
				smdict[mainkey] = subdict
			else:
				if smdict[mainkey][subkey]!=0.0:
					print 'error::line: %s' % line
					print 'error::subdict[%s] not zero: %f' % (subkey, subdict[subkey])
					return
				# before normalized
				smdict[mainkey][subkey] = rawcount # mainkey + subkey are unique in .2 and .1 files

	#process output
	smlist = [(k,sum(smdict[k].itervalues())) for k in smdict]
	totaltotal = sum([sum(smdict[k].itervalues()) for k in smdict])
	print 'total total: %.1f' % totaltotal
	smlist.sort(key=operator.itemgetter(1), reverse=True) # sort based on sum of nfreq to get rank i

	# init dist xtick (key) order
	sk = []
	print 'preparing output ..'
	for i in xrange(0, len(cm.aas01)):
		for j in xrange(0, len(cm.aas01)):
			sk.append('%s%s' % (cm.aas01[i], cm.aas01[j]))
	print repr(sk)
	print len(sk)

	#
	fout = open(outfile, 'w')
	for i in xrange(0, len(smlist)):
		mk, nfreq = smlist[i]
		h = sentropy(smdict[mk])

		dist = ' '.join([('%.6f' % (smdict[mk][k]/totaltotal)) for k in sk])
		outstr = '%d %s %.6f %.6f %s\n' % (i+1, mk, (nfreq/totaltotal), h, dist)
		fout.write(outstr)

	
	print 'save to %s.' % outfile
	fout.close()


# main entrance for pair substitution frequency stats
# dependence: proc_msa_pairsubstitution.py -> .pairfreq
def main():
	dispatch = {
		'pairfreqsum':pairfreqsum,
		'composition':composition,
		'ssentropy':ssentropy,
		'dist210':dist210,
		#'memo': proc_memo
	}

	if len(sys.argv)<2 or (sys.argv[1] not in dispatch):
		print 'Invalid cmd string'
		exit(1)
	else:
		dispatch[sys.argv[1]]()	
if __name__ == '__main__':
	main()