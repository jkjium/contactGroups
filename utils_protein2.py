import sys
from protein import protein
from AAmap import AAmap
import os.path
from collections import defaultdict
import numpy as np
import commp as cp

def cgfreq(arglist):
	if len(arglist) < 2:
		cp._err('Usage: python utils_protein2.py cgfreq cglistfile outfile.cgf')

	cglistfile = arglist[0]
	outfile = arglist[1]

	cgfreq = defaultdict(lambda:0)
	total = 0
	with open(cglistfile) as fp:
		for line in fp:
			# 2 R 26 D
			line = line.strip()
			if len(line) == 0:
				continue
			sarr = line.split(' ')
			key = '%s%s' % (sarr[1], sarr[3]) if sarr[1] <=sarr[3] else '%s%s' % (sarr[3], sarr[1])
			cgfreq[key]+=1
			total+=1
	'''
	print total, repr(cgfreq)
	tt=0
	for k in cgfreq:
		tt+= cgfreq[k]
	ttt=0
	'''
	cgf = [] # output with the same order
	for i in xrange(0, len(cp.aas01)):
		for j in xrange(i, len(cp.aas01)):
			key = '%s%s' % (cp.aas01[i], cp.aas01[j])
			cgf.append(float(cgfreq[key])/float(total))
			#ttt+=cgfreq[key]
			#print key,cgfreq[key],float(cgfreq[key])/float(total)
	#print total, tt, ttt

	with open(outfile, 'w') as fout:
		fout.write('%s\n' % (' '.join(['%.4f' % (f) for f in cgf])))
	cp._info('save to %s' % outfile)


def ncgfreq(arglist):
	if len(arglist)<3:
		cp._err('Usage: python utils_protein2.py ncgfreq cglistfile bgfreqfile outfile')

	cglistfile = arglist[0]
	bgfreqfile = arglist[1]
	outfile = arglist[2]

	cgfreq = defaultdict(lambda:0)
	with open(cglistfile) as fp:
		for line in fp:
			# 2 R 26 D
			line = line.strip()
			if len(line) == 0:
				continue
			sarr = line.split(' ')
			key = '%s%s' % (sarr[1], sarr[3]) if sarr[1] <=sarr[3] else '%s%s' % (sarr[3], sarr[1])
			cgfreq[key]+=1
	#print repr(cgfreq)

	bgfreq = np.loadtxt(bgfreqfile, delimiter=' ')
	#print repr(bgfreq)
	bgfreq = bgfreq/bgfreq.sum()
	#print repr(bgfreq)

	ncgfreq = []
	for i in xrange(0, len(cp.aas01)):
		for j in xrange(i, len(cp.aas01)):
			key = '%s%s' % (cp.aas01[i], cp.aas01[j])
			ncgfreq.append(cgfreq[key]/(bgfreq[i]*bgfreq[j]))
			#print key,cgfreq[key],ncgfreq
	nsum = sum(ncgfreq)
	#print nsum

	with open(outfile, 'w') as fout:
		fout.write(' '.join(['%.4f' % (f/nsum) for f in ncgfreq]))
	cp._info('save to %s' % outfile)

def seqaafreq(arglist):
	if len(arglist) < 2:
		cp._err('Usage: python utils_protein2.py seqaafreq 1t3r.pdb all (output: 1t3r.pdb.all.saf)')
	pdbfile = arglist[0]
	chainid = arglist[1]

	p = protein(pdbfile, chain=chainid)
	seqaa = p.seq.upper()
	aafreq = defaultdict(lambda:0)
	for aa in seqaa: 
		aafreq[aa]+=1

	outstr = ' '.join([('%d' % aafreq[aa]) for aa in cp.aas01])
	outfile = '%s.%s.saf' % (pdbfile, chainid)
	with open(outfile ,'w') as fout:
		fout.write('%s\n' % outstr)
	cp._info('save to %s' % outfile)

# write cutoff contact by specific method
# output: 1. method residue pdb file
# 		  2. contact file	
#def contactbycutoff(arglist):
def writecontact(arglist):
	if len(arglist) < 4:
		cp._err('Usage: python utils_protein2.py contactbycutoff 1t3r.pdb chain sgc cutoff')
	pdbfile = arglist[0] #sys.argv[2]
	chainid = arglist[1] #sys.argv[3]
	method = arglist[2] #sys.argv[4]
	cutoff = float(arglist[3]) #float(sys.argv[5])

	if method not in ['sgc', 'tip', 'ca']:
		cp._err('invalid method: %s' % method)

	p = protein(pdbfile, chain=chainid)
	if method =='sgc':
		ralist = p.atomsbyscgmcenter()
	elif method == 'tip':
		ralist = p.atomsbytip('AAtips.py')
	elif method == 'ca':
		ralist = p.ca
	# continue ...
	cglist = p.contactbycutoff(ralist, cutoff)

	outrafile = '%s.%s.%s' % (pdbfile, chainid, method)
	with open(outrafile, 'w') as fp:
		for a in ralist:
			fp.write(a.writeAtom())
	cp._info('save residue pdb: %s' % outrafile)

	outcgfile = '%s.%s.%s.cg' % (pdbfile, chainid, method)
	with open(outcgfile, 'w') as fp:
		for a,b in cglist:
			fp.write('%d %s %d %s\n' % (a.resSeq, cp.aa2a[a.resName], b.resSeq, cp.aa2a[b.resName]))
	cp._info('save cg file: %s' % outcgfile)

def foo(arglist):
	cp._info(repr(arglist))

if __name__ == '__main__':
	cp.dispatch(__name__)