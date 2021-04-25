import sys
from protein import protein
from AAmap import AAmap
import os.path
from collections import defaultdict
import numpy as np
import commp as cp

def alterbfactor(arglist):
	if len(arglist) < 4:
		cp._err('Usage:python utils_protein2.py alterbfactor pdbfile newbfactor.vec outfile ca_bf_or_not{1|0}')
	pdbfile = arglist[0]
	nbfile = arglist[1]
	outfile = arglist[2]
	cabfflag = int(arglist[3])

	rb = np.loadtxt(nbfile)
	nb = (rb - rb.min()) / (rb.max() - rb.min())

	p = protein(pdbfile)
	fout=open(outfile, 'w')

	if cabfflag == 1: # output allatom pdb
		if len(p.ca)!=len(nb):
			cp._err('length does not match. %d atoms, %d new bfactor values' % (len(p.ca), len(nb)))\
		# alter bfactor for all ca atoms
		for i in xrange(0, len(p.ca)):
			a=p.ca[i]
			a.tempFactor=nb[i]
		# same all atom pdb
		for a in p.atoms:
			if 'CA' not in a.name:
				a.tempFactor=0
			fout.write(a.writeAtom())
	else:
		if len(p.atoms)!=len(nb):
			cp._err('length does not match. %d atoms, %d new bfactor values' % (len(p.atoms), len(nb)))
		# alter bfactor
		for i in xrange(0, len(p.atoms)):
			a= p.atoms[i]
			a.tempFactor = nb[i]
			fout.write(a.writeAtom())		

	fout.close()
	cp._info('save new pdb to %s' % outfile)

# alter b factor according to the given resi
def alterbresi(args):
	assert len(args) == 4, 'Usage: python utils_protein2.py alterbresi pdbfile newbfactor.vec resifile.tick outfile'
	pdbfile = args[0]
	nbfile = args[1]
	resifile = args[2]
	outfile = args[3]

	rb = np.loadtxt(nbfile)
	nb = (rb - rb.min()) / (rb.max() - rb.min()) # normalize to 0 - 1

	resilist = cp.loadlines(resifile)
	assert len(resilist) == len(nb)

	nbdict = dict((int(resilist[i]), nb[i]) for i in range(len(resilist)))

	p = protein(pdbfile)
	fout=open(outfile, 'w')

	# alter bfactor for all ca atoms
	for i in xrange(0, len(p.ca)):
		a=p.ca[i]
		a.tempFactor=nbdict[a.resSeq] if a.resSeq in nbdict else 0
	# same all atom pdb
	for a in p.atoms:
		if 'CA' not in a.name:
			a.tempFactor=0
		fout.write(a.writeAtom())
	fout.close()
	cp._info('save new pdb to %s' % outfile)



# alter b factor according to the given (resi, b) .vec2 file
# b factor normalization should be done by upper stream processing
def alterbresi2(args):
	assert len(args) == 3, 'Usage: python utils_protein2.py alterbresi2 pdbfile resi_nb.vec2 outpdbfile'
	pdbfile = args[0]
	nbfile = args[1]
	outfile = args[2]

	resi_b = lambda x: (int(x[0]), float(x[1]))
	nbdict = dict(resi_b(line.split()) for line in cp.loadlines(nbfile))

	p = protein(pdbfile)
	fout=open(outfile, 'w')

	# alter bfactor for all ca atoms
	for i in xrange(0, len(p.ca)):
		a=p.ca[i]
		a.tempFactor=nbdict[a.resSeq] if a.resSeq in nbdict else 0
	# same all atom pdb
	for a in p.atoms:
		if 'CA' not in a.name:
			a.tempFactor=0
		fout.write(a.writeAtom())
	fout.close()
	cp._info('save new pdb to %s' % outfile)	


# replace resSeq from oldresi to new resi
# nresifile a single column data contains the new resi
# oresifile gives which old resi to be replaced
def alterresseq(args):
	assert len(args) == 4, 'Usage: python utils_protein2.py alterresseq pdbfile newresi.vec oldresi.vec(tick) outpdbfile'
	pdbfile = args[0]
	nresifile =args[1]
	oresifile =args[2]
	outfile = args[3]

	nresi = np.loadtxt(nresifile)
	oresi = np.loadtxt(oresifile)
	assert len(nresi) == len(oresi), 'Error, length mismatch: nresi: %d, oresi: %d' % (len(nresi), len(oresi))

	nrdict = dict((int(oresi[i]), int(nresi[i])) for i in range(len(nresi)))
	count=0
	p = protein(pdbfile)
	with open(outfile, 'w') as fout:
		for a in p.atoms:
			if a.resSeq in nrdict:
				a.resSeq = nrdict[a.resSeq] 
				if 'CA' in a.name:
					count+=1
			fout.write(a.writeAtom())
	cp._info('save new resi to %d residues in new pdb file %s' % (count, outfile))


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

	cgf = [] # output with the same order
	for i in xrange(0, len(cp.aas01)):
		for j in xrange(i, len(cp.aas01)):
			key = '%s%s' % (cp.aas01[i], cp.aas01[j])
			cgf.append(float(cgfreq[key])/float(total))

	with open(outfile, 'w') as fout:
		fout.write('%s\n' % (' '.join(['%.4f' % (f) for f in cgf])))
	cp._info('save to %s' % outfile)


# output reisdue contact pairs for two input pdbs
# input: pdb1 pdb2
# output: pdb1 {chain resi resn} pdb2 {chain resi resn} mindist_of_all_atoms
def pdbcontact(args):
	assert len(args)==3, 'Usage: python utils_protein2.py pdbcontact pdbfile1 pdbfile2 outfile'
	p1 = protein(args[0])
	p2 = protein(args[1])
	outfile = args[2]
	outlist = []
	for r1 in p1.resAtomsDict:
		for r2 in p2.resAtomsDict:
			res1atom = p1.resAtomsDict[r1][0]
			res2atom = p2.resAtomsDict[r2][0]
			mindist = min([np.linalg.norm(np.array((a.x, a.y, a.z))-np.array((b.x, b.y, b.z))) for a in p1.resAtomsDict[r1] for b in p2.resAtomsDict[r2]])
			outstr = '%s %d %s %s %d %s %.4f' % (res1atom.chainID, res1atom.resSeq, cp.aa2a[res1atom.resName], res2atom.chainID, res2atom.resSeq, cp.aa2a[res2atom.resName], mindist)
			outlist.append(outstr)
	with open(outfile, 'w') as fout:
		fout.write('%s\n' % '\n'.join(outlist))
	cp._info('save contacts to %s' % outfile)


# input: pdb file, target residue IDs, cutoff distance
# output: a flat file
# format:
# resi n1,n2,n3... neighborstring
# atom.self.resSeq (int)
# protein.contactbyallatom(self, chain, resi, cutoff, seqcutoff=0.0); return list: ['A161','A161', A256', ...]
# example:
# $ python utils_protein2.py neighborsflatline 1gzh_A.pdb A161,A255 4.0
# A161 A255 FVRAMIVIYLTITRNF A109,A157,A158,A159,A160,A162,A173,A195,A234,A252,A253,A254,A256,A267,A268,A270
def neighborsflatline(arglist):
	if len(arglist) < 3:
		cp._err('Usage: python utils_protein2.py neighborsflatline pdbfile A2,B13 4.0')

	pdbfile = arglist[0]
	p = protein(pdbfile)
	reslist = [(s[0], s[1:]) for s in arglist[1].split(',')]
	targetstr = arglist[1].replace(',',' ')
	cutoff = float(arglist[2])

	# get neighbor information
	n=[]
	# [('A', '161'),('A','161'), ('A','256'), ...]
	for r in reslist:
		n+=p.contactbyallatom(r[0], int(r[1]), cutoff)
	# remove redundancy
	neighbors = list(set(n))
	neighbors.sort()
	resstr = ','.join(neighbors)
	# get local neighbors sequence
	# p.resDict = 'B529': (132, 'V')
	neighborseq = ''.join([p.resDict[resi][1] for resi in neighbors])
	print '%s %s %s' % (targetstr, neighborseq, resstr)

def writechain(arglist):
	if len(arglist) < 3:
		cp._err('Usage: python utils_protein2.py writechain pdbfile chainID outfile')
	pdbfile = arglist[0]
	c = arglist[1]
	outfile = arglist[2]

	p = protein(pdbfile, chain=c)

	with open(outfile, 'w') as fout:
		for a in p.atoms:
			fout.write(a.writeAtom())
	cp._info('save %s chain %s to %s' % (pdbfile, c, outfile))



# input: pdb file, target residue IDs, cutoff distance
# output: a flat file
# format:
# resi n1,n2,n3... neighborstring
# atom.self.resSeq (int)
# protein.contactbyallatom(self, chain, resi, cutoff, seqcutoff=0.0); return list: ['A161','A161', A256', ...]
# example:
# $ python utils_protein2.py neighborsflatline 1gzh_A.pdb A161,A255 4.0
# A161 A255 FVRAMIVIYLTITRNF A109,A157,A158,A159,A160,A162,A173,A195,A234,A252,A253,A254,A256,A267,A268,A270
def neighborsflatline2(arglist):
	if len(arglist) < 4:
		cp._err('Usage: python utils_protein2.py neighborsflatline pdbfile A2,B13 opt{all|ca} 4.0')

	pdbfile = arglist[0]
	p = protein(pdbfile)
	reslist = [(s[0], s[1:]) for s in arglist[1].split(',')]
	targetstr = arglist[1].replace(',',' ')
	opt = arglist[2]
	cutoff = float(arglist[3])

	# get neighbor information
	n=[]
	# [('A', '161'),('A','161'), ('A','256'), ...]
	for r in reslist:
		if opt=='ca':
			n+=p.neighborsbyca(r[0], int(r[1]), cutoff)
		else:
			n+=p.contactbyallatom(r[0], int(r[1]), cutoff)
	# remove redundancy
	neighbors = list(set(n))
	neighbors.sort()
	resstr = ','.join(neighbors)
	# get local neighbors sequence
	# p.resDict = 'B529': (132, 'V')
	neighborseq = ''.join([p.resDict[resi][1] for resi in neighbors])
	print '%s %s %s' % (targetstr, neighborseq, resstr)


# input .pdb file
# output: a flat file, recording the minimum distance among all atoms of the current pair of residues
# A123 A124 5.356
def writeresdists(args):
	# residistbyallatom
	assert len(args) == 2, 'Usage: python utils_protein2.py writeresdists pdbfile outfile'
	infile = args[0]
	outfile = args[1]

	p = protein(infile)
	# [(rA, rB, dist), (.), ..]
	# A < B
	dtuplelist = p.residistbyallatom()
	with open(outfile, 'w') as fout:
		fout.write('%s\n' % ('\n'.join(['%s %s %s %s %.4f' % (t[0][0], t[0][1:], t[1][0], t[1][1:], t[2]) for t in dtuplelist ]))) 
	cp._info('save to %s' % outfile)



def ncgfreq(arglist):
	if len(arglist)<3:
		cp._err('Usage: python utils_protein2.py ncgfreq cglistfile bgfreqfile outprefix')

	cglistfile = arglist[0]
	bgfreqfile = arglist[1] # list of all bg freq
	outfile = arglist[2]

	# count all the cg observations
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

	# sum up AA count 
	aafreq = np.loadtxt(bgfreqfile, delimiter=' ')
	bgfreq = np.sum(aafreq, axis=0)

	bgsum = float(sum(bgfreq))
	# save .allbgfreq
	with open(outfile+'.bgfreq', 'w') as fout:
		fout.write(' '.join(['%.4f' % (i/bgsum) for i in bgfreq]))
	cp._info('save background AA frequency to %s (sum %d)' % (outfile+'.bgfreq', bgsum ))

	bgfreq = bgfreq/bgsum
	ncgfreq = [] # normalized cg freq
	ocgfreq = [] # orginal cg freq
	for i in xrange(0, len(cp.aas01)):
		for j in xrange(i, len(cp.aas01)):
			key = '%s%s' % (cp.aas01[i], cp.aas01[j])
			ncgfreq.append(cgfreq[key]/(bgfreq[i]*bgfreq[j]))
			ocgfreq.append(cgfreq[key])

	osum = float(sum(ocgfreq))
	with open(outfile+'.ocgfreq', 'w') as fout:
		fout.write(' '.join(['%.4f' % (f/osum) for f in ocgfreq]))
	cp._info('save original cg frequency to %s' % (outfile+'.ocgfreq'))

	nsum = float(sum(ncgfreq))
	with open(outfile+'.ncgfreq', 'w') as fout:
		fout.write(' '.join(['%.4f' % (f/nsum) for f in ncgfreq]))
	cp._info('save to normalized cg frequency %s' % (outfile+'.ncgfreq'))


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


# cut pdb using one sequence segment
def splitpdbbyseq(arglist):
	if len(arglist) < 3:
		cp._err('Usage: python utils_protein2.py splitpdbbyseq pdbfile seqfafile outpdbfile')
	pdbfile = arglist[0]
	seqfile = arglist[1]
	outfile = arglist[2]

	# load pdb
	p = protein(pdbfile)
	# load seqfa
	for h,s in cp.fasta_iter(seqfile):
		seqheader = h
		seqbody=s.translate(None, ''.join(cp.gaps)).upper()
	
	outres = p.slicebyseq(seqbody)
	if len(outres) == 0:
		cp._err('mismatch %s %s' % (pdbfile, seqfile))
	with open(outfile, 'w') as fout:
		for r in outres:
			for a in r:
				fout.write(a.writeAtom())
	cp._info('save %d residues to %s' % (len(seqbody), outfile))


# write cutoff contact by specific method
# output: 1. method residue pdb file
# 		  2. contact file	
#def contactbycutoff(arglist):
def writecontact(arglist):
	if len(arglist) < 4:
		cp._err('Usage: python utils_protein2.py writecontact 1t3r.pdb chain sgc cutoff')
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
			fp.write('%d %s %d %s chain %s %s\n' % (a.resSeq, cp.aa2a[a.resName], b.resSeq, cp.aa2a[b.resName], a.chainID, b.chainID))
	cp._info('save cg file: %s' % outcgfile)


# input: atom list
# output: pdist between every two atoms
# format: {(chainA, resiA), (chainB, resiB), dist}
def _respdist(atomset):
	pdistlist = []
	for i in range(0, len(atomset)):
		for j in range(i+1, len(atomset)):
			a = atomset[i]
			b = atomset[j]
			chainA = a.chainID
			chainB = b.chainID
			rA = a.resSeq
			rB = b.resSeq
			dist =  np.linalg.norm(np.array((a.x, a.y, a.z))-np.array((b.x, b.y, b.z)))
			pdistlist.append(('%s%d' % (chainA,rA), '%s%d' % (chainB, rB), dist))
	return pdistlist


# input protein instance, pdist evaluation method {sgc, tip, ca, all}
# return {chain+res1,chain+res2,pdist}
def _resdist_by_method(p, method):
	cp._info('calculating pdist using %s atoms' % method)
	if method =='sgc':
		dtuplelist = _respdist(p.atomsbyscgmcenter())
	elif method == 'tip':
		dtuplelist = _respdist(p.atomsbytip())
	elif method == 'ca':
		dtuplelist = _respdist(p.ca)
	else:
		# [(chainA, rA, chainB, rB, dist), (.), ..]
		dtuplelist = p.residistbyallatom()
	print dtuplelist[0]
	return dtuplelist

# write cutoff contact by specific method
# output: 1. method residue pdb file
# 		  2. contact file	
#def contactbycutoff(arglist):
def writedistby(args):
	assert len(args) == 3, 'Usage: python utils_protein2.py writedistby pdbfile method={tip, sgc, tip, all} outfile'
	pdbfile = args[0]
	method = args[1]
	outfile = args[2]

	dtuplelist = _resdist_by_method(protein(pdbfile), method)
	with open(outfile, 'w') as fout:
		fout.write('%s\n' % ('\n'.join(['%s %s %s %s %.4f' % (t[0][0], t[0][1:], t[1][0], t[1][1:], t[2]) for t in dtuplelist ])))
	cp._info('save %d pdist to %s' % (len(dtuplelist),outfile))


# for protsub
# write all chain sequence into one outfile
def writefafile(arglist):
	if len(arglist) < 2:
		cp._err('Usage: python utils_protein2.py writefafile pdbfile outfile')
	pdbfile = arglist[0]
	outfile = arglist[1]

	p = protein(pdbfile, chain='all')
	fout = open(outfile, 'w')
	fout.write('>%s/0-%d\n%s\n' % (outfile, len(p.seq), p.seq.lower()))
	fout.close()
	cp._info('write sequence to %s, len %d' % (outfile, len(p.seq)))



def writeseqfa(arglist):
	if len(arglist) == 2:
		chainid = arglist[1]
		outfile = '%s.%s.fa' % (arglist[0], chainid)
		head = '%s.%s' % (arglist[0], chainid)
	elif len(arglist) == 1:
		chainid = 'all'
		outfile = '%s.fa' % arglist[0]
		head = arglist[0]
	else:
		cp._err('Usage: python utils_protein.py writeseqfa 1t3r.pdb {A}')

	pdbfile = arglist[0]

	p = protein(pdbfile, chain=chainid)
	fout = open(outfile, 'w')
	fout.write('>%s/0-%d\n%s\n' % (head, len(p.seq), p.seq.lower()))
	fout.close()
	cp._info('write sequence to %s, len %d' % (outfile, len(p.seq)))

def writeseq(arglist):
	if len(arglist) == 2:
		chainid = arglist[1]
		outfile = '%s.%s.seq' % (arglist[0], chainid)
		head = '%s.%s' % (arglist[0], chainid)
	elif len(arglist) == 1:
		chainid = 'all'
		outfile = '%s.seq' % arglist[0]
		head = arglist[0]
	else:
		cp._err('Usage: python utils_protein.py writeseqfa 1t3r.pdb {A}')

	pdbfile = arglist[0]

	p = protein(pdbfile, chain=chainid)
	fout = open(outfile, 'w')
	fout.write('%s\n' % (p.seq.lower()))
	fout.close()
	cp._info('write sequence to %s, len %d' % (outfile, len(p.seq)))

def foo(arglist):
	cp._info(repr(arglist))

if __name__ == '__main__':
	cp.dispatch(__name__)
