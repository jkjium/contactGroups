'''
get resid list of varname
'''
import sys
from protein import protein
from AAmap import AAmap
import commp as cp
import os.path
import numpy as np

def resn2bfactor():
	if len(sys.argv) < 3:
		print 'resn2bfactor(): replace b factor values with residue type.'
		print 'resn2bfactor(): used for pymol spectrum b'
		return
	scoreValue = {
							'X':0,'-': 0,'.': 0,'A': 1,'C': 2,'D': 3,'E': 4,'F': 5,'G': 6,'H': 7,'I': 8,'K': 9,
							'L': 10,'M': 11,'N': 12,'P': 13,'Q': 14,'R': 15,'S': 16,'T': 17,'V': 18,'W': 19,'Y': 20, 'B': 3
						}
	aamap = AAmap()

	pdbfile = sys.argv[2]
	p = protein(pdbfile)
	outfile = '%s_rb.pdb' % pdbfile[:-4]
	fout = open(outfile, 'w')
	for a in p.atoms:
		newBFactor = scoreValue[aamap.getAAmap(a.resName)]
		print 'new b-factor: [%s : %s] -> %d' % (a.resName, aamap.getAAmap(a.resName), newBFactor)
		a.tempFactor = newBFactor
		fout.write(a.writeAtom())
	fout.close()
	print 'Output file: %s' % outfile


def pdbcut():
	if len(sys.argv) < 5:
		print 'pdbcut(): write pdb by residue segment'
		print 'pdbcut(): python utils_protein.py pdbcut 1t3r.pdb A 5-15'
		print 'pdbcut(): python utils_protein.py pdbcut 1t3r.pdb all 5-15'
		return

	pdbfile = sys.argv[2]
	chain = sys.argv[3]
	rangeStr = sys.argv[4]


	pdbname = pdbfile[0:4]

	print 'pdbcut():pdbfile: %s' % pdbfile
	print 'pdbcut():pdb: %s' % pdbname
	print 'pdbcut():chain: %s' % chain

	p = protein(pdbfile)
	if rangeStr != 'all':
		rangeArray = rangeStr.split('-')
		rBegin = int(rangeArray[0])
		rEnd = int(rangeArray[1])
		outfile = '%s-%s-%d-%d.rpdb' % (pdbname, chain, rBegin, rEnd)
	else:
		rBegin = -999
		rEnd = 999999
		outfile = '%s-%s.rpdb' % (pdbname, chain)
	print 'pdbcut():residue range: %s, %d - %d' % (rangeStr, rBegin, rEnd)

	out = []
	if chain == 'all':
		for a in p.atoms:
			if a.resSeq <= rEnd and a.resSeq >= rBegin:
				out.append(a)
	else:
		for a in p.atoms:
			if (a.resSeq <= rEnd and a.resSeq >= rBegin and a.chainID.lower() == chain.lower()):
				out.append(a)

	fout = open(outfile, 'w')
	print 'pdbcut():output: %s' % outfile
	print 'pdbcut():%s: %d atoms written.' % (outfile, len(out))
	for a in out:
		fout.write(a.writeAtom())
	fout.close()


# check sgc, tip, ca atoms
def pdbscreen():
	if len(sys.argv) < 4:
		cp._err('Usage: python utils_protein.py pdbscreen pdbfile all|{A}')

	pdbfile = sys.argv[2]
	chainid = sys.argv[3]

	if os.path.isfile(pdbfile):
		p = protein(pdbfile, chain=chainid)
		#report = np.array([1.0*len(p.atomsbyscgmcenter())/len(p.resDict), 1.0*len(p.atomsbytip())/len(p.resDict), 1.0*len(p.ca)/len(p.resDict)])
		report = np.array([1.0*len(p.atomsbyscgmcenter())/len(p.seq), 1.0*len(p.atomsbytip())/len(p.seq), 1.0*len(p.ca)/len(p.seq)])

		status = 0
		# no information available
		if report.sum()== 0:
			status = -2
		elif (report[0] == report[1]) and (report[1] == report[2]):
			status = 0
		else:
			# information partially available
			status = 1

		cp._info('%s %s stat %d %d %.3f %.3f %.3f' % (pdbfile, chainid, status, len(p.seq), report[0], report[1], report[2]))
	else:
		# file does not exist
		cp._info('%s %s stat -1 0 0 0 0' % (pdbfile, chainid))


def writeseq():
	if len(sys.argv) == 4:
		chainid = sys.argv[3]
		outfile = '%s.%s.seq' % (sys.argv[2], chainid)
	elif len(sys.argv) == 3:
		chainid = 'all'
		outfile = '%s.seq' % sys.argv[2]
	else:
		print 'Usage: python utils_protein.py writeseq 1t3r.pdb {A}'
		print 'output: 1t3r.pdb.{A.}seq'
		return

	pdbfile = sys.argv[2]


	p = protein(pdbfile, chain=chainid)
	fout = open(outfile, 'w')
	fout.write(p.seq.lower()+'\n')
	fout.close()
	print 'writeseq(): outfile: %s' % outfile


def writeseqfa():
	if len(sys.argv) == 4:
		chainid = sys.argv[3]
		outfile = '%s.%s.fa' % (sys.argv[2], chainid)
		head = '%s.%s' % (sys.argv[2], chainid)
	elif len(sys.argv) == 3:
		chainid = 'all'
		outfile = '%s.fa' % sys.argv[2]
		head = sys.argv[2]
	else:
		print 'Usage: python utils_protein.py writeseq 1t3r.pdb {A}'
		print 'output: 1t3r.pdb.{A.}seq'
		return

	pdbfile = sys.argv[2]


	p = protein(pdbfile, chain=chainid)
	fout = open(outfile, 'w')
	fout.write('>%s\n%s\n' % (head, p.seq.lower()))
	fout.close()
	print 'writeseq(): outfile: %s' % outfile


# irite tip pdb
def writetip():
	if len(sys.argv) < 3:
		print 'writetip(): write tip pdb file'
		print 'writeseq(): python utils_protein.py writetip pdbfile'
		print 'writetip(): output: 1t3r.tip'
		return

	pdbfile = sys.argv[2]
	outfile = pdbfile+'.tip'

	p = protein(pdbfile)
	p.writeTips('AAtips.py',outfile)
	print 'save to %s' % outfile

# side chain geom center
def writesgc():
	if len(sys.argv) < 3:
		print 'writeseq(): python utils_protein.py writesgc pdbfile'
		print 'writetip(): output: 1t3r.pdb.sgc'
		return

	pdbfile = sys.argv[2]
	outfile = sys.argv[2] + '.sgc'

	fout = open(outfile, 'w')
	p = protein(pdbfile)
	for a in p.atomsbyscgmcenter():
		fout.write(a.writeAtom())
	fout.close()
	print 'save to %s' % outfile


def writeca():
	if len(sys.argv) < 3:
		print 'writeca(): python utils_protein.py writeca pdbfile {chain}'
		print 'writeca(): output: 1k2p.pdb.chain.ca'
		return
	chain = 'all'
	pdbfile = sys.argv[2]
	if len(sys.argv) == 4:
		chain = sys.argv[3]

	outfile = pdbfile + '.ca' if chain == 'all' else '%s.%s.ca' % (pdbfile, chain)
	p = protein(pdbfile)
	p.writeCA(outfile, chain)
	print 'save to %s' % outfile


def dumpseqflat():
	if len(sys.argv) < 4:
		print 'Printout sequence in flat foramt: [seq name] [sequence]'
		print 'python utils_protein.py dumpseqflat pdbfile chain'
		print 'python utils_protein.py dumpseqflat 2gag.pdb A'
		return

	pdbfile = sys.argv[2]
	chain = sys.argv[3]
	p = protein(pdbfile, chain=chain)
	print '%d %s %s' % (len(p.seq), pdbfile, p.seq)

# write cutoff contact by specific method
# output: 1. method residue pdb file
# 		  2. contact file	
def contactbycutoff():
	if len(sys.argv) < 6:
		cp._err('Usage: python utils_protein.py contactbycutoff 1t3r.pdb chain sgc cutoff')
	pdbfile = sys.argv[2]
	chainid = sys.argv[3]
	method = sys.argv[4]
	cutoff = float(sys.argv[5])

	if method not in ['sgc', 'tip', 'ca']:
		cp._err('invalid method: %s' % method)

	p = protein(pdbfile, chain=chainid)
	if method =='sgc':
		ralist = p.atomsbyscgmcenter()
	elif method == 'tip':
		ralist = p.atomsbytip('AAtips.py')
	# elif method == 'ca':
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


def main():
	if len(sys.argv)<2:
		print 'Usage: python utils_protein.py cmd pdbfile [args ...]'
		return

	dispatch = {
		'resn2bfactor': resn2bfactor, 'pdbcut': pdbcut, 'writeseq':writeseq, 'writetip':writetip, 'dumpseqflat':dumpseqflat,
		'writeca':writeca, 'writesgc':writesgc, 'writeseqfa':writeseqfa,
		'contactbycutoff':contactbycutoff,
		'pdbscreen': pdbscreen
	}

	cmd = sys.argv[1]

	if sys.argv[1] not in dispatch:
		print 'invalid cmd: %s' % sys.argv[1]
		return

	for key in dispatch:
		if key == cmd:
			dispatch[key]()


if __name__ == '__main__':
	main()
