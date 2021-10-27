import commp as cp
import numpy as np

# input: 
#	1. filtered spike protein sequences
# 	2. stub file contains timestamp needed
# output: sampled sequences (for later dendrogram)
def timesample(args):
	infile = args[0]
	stubfile = args[1]
	outfile = args[2]

	stub = cp.loadlines(stubfile)
	n = len(stub)
	
	outdict ={}
	for s in cp.fasta_iter(infile):
		sarr = s[0].split('|') 
		if sarr[2] in stub:
			outdict[sarr[2]] = '>%s\n%s' % (s[0], s[1])
			stub.remove(sarr[2])
		if len(stub) == 0:
			break
	cp._info('%d/%d sequences sampled.' % (len(outdict), n))
	ts = outdict.keys()
	ts.sort()
	fout= open(outfile, 'w')
	with open(outfile, 'w') as fout:
		fout.write('%s\n' %  ('\n'.join([outdict[fa] for fa in ts])))
	cp._info('save sampled fas to %s' % outfile)



# screen 600,000 sprotein_nr sequences for distributions of length and freq(X's)
# input:
# bbmb_kjia_admin@BB-RJER-3630519 ~/workspace/covid-19/all/sprotein 2021-10-25 11:37:24
# sprotein_nr.fasta
# output: .vec2 {1.length 2.n.Xs}
def seqscreen(args):
	infile = args[0]
	outfile = args[0]+'.stat'
	fout = open(outfile, 'w')
	count=0
	for s in cp.fasta_iter(infile):
		xc=0
		for i in s[1]:
			if 'X' == i:
				xc+=1
		fout.write('%d %d\n' % (len(s[1]), xc))
		count+=1
		if count%10000==0:
			print('%d' % count)
	fout.close()
	cp._info('save to %s' % outfile)

# filter sequences by length cutoff and freq{X}
def seqfilter(args):
	infile =args[0]
	lencutoff = int(args[1])
	xcutoff = int(args[2])
	outfile = '%s.%d.fa' % (infile, lencutoff) 
	fout=open(outfile, 'w')
	count=0
	for s in cp.fasta_iter(infile):
		xc=0
		for i in s[1]:
			if 'X' == i:
				xc+=1		
		if len(s[1])<lencutoff or xc>=xcutoff:
		#if 'X' in s[1] and len(s[1])<lencutoff:
			continue
		fout.write('>%s\n%s\n' % (s[0], s[1]))
		count+=1
		if count%10000==0:
			print(count)
	fout.close()
	cp._info('save %d sequences to %s' % (count, outfile))


# return an np matrix n (number of positions) x 20 (20 types of amino acid)
def logomat(args):
	assert len(args) == 2, 'Usage: python proc_covid.py logomat in.tsv out.mat'
	from collections import defaultdict
	infile = args[0]
	outprefix = args[1]
	posdict = {}
	entdict = defaultdict(list)
	for s in cp.loadtuples(infile):
		p = int(s[1])
		entdict[p].append(cp.aascore['aa'][s[0]])
		entdict[p].append(cp.aascore['aa'][s[2]])
		#print(s)
		if(p not in posdict):
			posdict[p] = defaultdict(int)
			posdict[p][s[0]]+=1
			posdict[p][s[2]]+=1
		else:
			posdict[p][s[0]]+=1
			posdict[p][s[2]]+=1

	klist = posdict.keys()
	klist.sort()
	#print(klist)
	outmat = []
	outent = []
	for i in klist:
		aafreqlist = [posdict[i][c] for c in cp.aas01]
		n= np.array(entdict[i])
		outent.append('%d %.4f' % (i, cp.entropy([n])))
		
		'''
		if(i==270):
			print aafreqlist
		'''
		outmat.append(aafreqlist)

	#print('%s\n' % (' '.join(cp.aas01)))
	#print(np.array(outmat))
	npoutmat = np.array(outmat)
	npoutcol = np.array(klist)

	outmatfile = '%s.mat' % outprefix
	np.savetxt(outmatfile, npoutmat, fmt='%d')
	outcolfile = '%s.col' % outprefix
	np.savetxt(outcolfile, npoutcol, fmt='%d')
	outentfile = '%s.ent' % outprefix
	with open(outentfile, 'w') as fout:
		fout.write('%s\n' % ('\n'.join(outent)))

	cp._info('save to %s %s %s' % (outmatfile, outcolfile, outentfile))

# input 1:  seqlogo.mat from logomat()
# input 2: col_entropy.vec2 file 'seqlogo.ent'
# output sub-mat, positon_list
def maxent(args):
	assert len(args) == 3, 'Usage:python proc_covid.py maxentropy inmatfile inentfile 30{segment length}'
	inmatfile = args[0]
	inentfile = args[1]
	s = int(args[2])
	logomat = np.loadtxt(inmatfile)
	entmat = np.loadtxt(inentfile)

	entlist = entmat[:,1]
	poslist = entmat[:,0]

	n = len(poslist)
	maxent = 0.
	print(n, s)
	for i in range(0,n,s):
		sument = np.sum(entlist[i:i+30])
		print(sument, maxent)
		if(sument>maxent):
			maxent = sument
			maxs = poslist[i:i+30]
		print(i,i+30)
		#print(maxs)
		print("")

	print('maxent: %.4f' % maxent)
	print('maxseg: %s' % (','.join(['%d' % i for i in maxs])))
	


	

# input: fasta file contains all the protein sequences of all the variants
def uniqvariantseq(args):
	assert len(args)==2, 'Usage: python proc_covid.py variants.fa out.fa'

	namemap = {'B.1.1.7': 'Alpha', 'B.1.351': 'Beta', 'B.1.617.2': 'Delta', 'AY.1': 'Delta', 
				'AY.2': 'Delta', 'AY.3': 'Delta', 'AY.3.1': 'Delta', 'P.1': 'Gamma', 
				'B.1.427': 'Epsilon', 'B.1.429': 'Epsilon', 'B.1.525': 'Eta', 'B.1.526': 'Iota', 
				'B.1.617.1': 'Kappa', 'B.1.617.3': 'na', 'sars2': 'na'} 

	infafile = args[0]
	outfile = args[1]

	seqdict = {}
	for h, s in cp.fasta_iter(infafile):
		if 'polyprotein' in h:
			cp._info('skipping %s' % h)
			continue
		sarr = h.split('|') # surface glycoprotein|B.1.525
		seqtype = sarr[0]
		seqname = sarr[1]
		varname = namemap[sarr[1]]
		if(s in seqdict):
			seqdict[s] = '%s|%s %s' % (seqdict[s], varname, seqname)
		else:
			seqdict[s] = '%s|%s %s' % (seqtype, varname, seqname)

	outstr = '\n'.join(['>%s\n%s' % (seqdict[s], s) for s in seqdict])
	with open(outfile, 'w') as fout:
		fout.write('%s\n' % outstr)
	cp._info('save fasta to %s.' % outfile)



def fa2flat(args):
	assert len(args)==2, 'Usage: python proc_fa2flat.py fa2flat in.fa outfile'
	msafile = args[0]
	outfile = args[1]
	
	fout=open(outfile, 'w')
	for head, seq in cp.fasta_iter(msafile):
		fout.write('%s %s\n' % (head.ljust(40), seq))
					
	fout.close()


if __name__=='__main__':
	cp.dispatch(__name__)