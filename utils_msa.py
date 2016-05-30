'''
get msa position id 
'''
import sys
import numpy as np
from msa import msa
from protein import protein

# given a residue number output the corresponding position in msa file
def resi2msai():
	if len(sys.argv) < 5:
		print 'resi2msai: given a residue number output the corresponding position in msa'
		print 'example:python utils_msa.py resi2msai 1k2p_PF07714_full.fa 1k2p.pdb B641\n'
		return

	msafile = sys.argv[2]
	pdbfile = sys.argv[3]
	resi = sys.argv[4]

	#print repr((msafile, pdbfile, resi))
	m = msa(msafile)
	p = protein(pdbfile)

	resIdx = p.resDict[resi]
	posMap = m.getPosMap(p)[0]

	for i in posMap:
		if i == resIdx[0]:
			print '[Res: %s] : (seqi: %d (%s) - msai: %d (%s))' % (resi, resIdx[0], resIdx[1], posMap[i], m.msaArray[0][1][posMap[i]])
			break


def sdii2resi():
	if len(sys.argv) < 5:
		print 'resi2target: given a residue number output the corresponding position in target msa'
		print 'example:python utils_msa.py sdii2resi PF07714_full.fa.r50 BTK_HUMAN 1k2p.pdb PF07714_full.fa.r50.3128_3_sdii\n'
		print 'output: PF07714_full.fa.r50.3128_3_sdii_resi'
		return

	msafile = sys.argv[2]
	target = sys.argv[3]
	pdbfile = sys.argv[4]
	sdiifile = sys.argv[5]

	print 'msafile: %s\ntarget header: %s\npdbfile: %s\nsdii file: %s' % (msafile, target, pdbfile, sdiifile)
	m = msa(msafile)
	p = protein(pdbfile)
	rtmap = m.getResiTargetMap(p, target)
	if len(rtmap) < 1:
		print 'error occoured in generating rtmap'
		return
	#print '%s: %s' % (tvar, repr(rtmap[tvar]))
	# construct trmap from rtmap
	# 3128: (B641, 'R')
	trmap = {}
	#trmap = {v: k for k, v in rtmap.iteritems()}
	for k in rtmap:
		msai, resn = rtmap[k]
		if msai in trmap:
			print 'error. duplicate key [%d] in rtmap' % msai
			return
		trmap[msai] = (k, resn)

	#print trmap

	# read sdii file
	with open(sdiifile) as f:
		sdiilines = f.readlines()

	outfile = sdiifile + '_resi'
	fout = open(outfile, 'w')

	# 52 [pid:20029] 926-3089-3128 0.001106226720675
	count = 0
	for line in sdiilines:
		count += 1
		print '%d/%d processed ...' % (count, len(sdiilines))
		strArr = line.strip().split(' ')
		msailist = strArr[2].split('-')
		sdiivalue = strArr[3]
		fout.write('%s %s\n' % ('-'.join([repr(trmap[int(i)]) for i in msailist]), sdiivalue))
	fout.close()
	print 'done.\noutput file: [%s]' % outfile

def msai2resi():
	pass

# output: B641: (3128, 'R')
def resi2target():
	if len(sys.argv) < 5:
		print 'resi2target: given a residue number output the corresponding position in target msa'
		print 'example:python utils_msa.py resi2target PF07714_full.fa.r50 BTK_HUMAN 1k2p.pdb B641\n'
		return

	msafile = sys.argv[2]
	target = sys.argv[3]
	pdbfile = sys.argv[4]
	tvar = sys.argv[5]

	print 'msafile: %s\ntarget header: %s\npdbfile: %s\ntarget variable: %s' % (msafile, target, pdbfile, tvar)
	m = msa(msafile)
	p = protein(pdbfile)
	print p.resDict[tvar]
	rtmap = m.getResiTargetMap(p, target)
	if len(rtmap) < 1:
		return
	print '%s: %s' % (tvar, repr(rtmap[tvar]))



# function for parsing sdii result
# msai -> seqi -> 'B529(V)'
# pdbseqDict: 132 : 'B529(V)'
# sdiline: [1042-2032-3128 0.006242240179705]
# discard for implant alignment
'''
def sdiiparse(sdiiline, msai2seqi, pdbseqDict):
	split1 = sdiiline.split(' ')
	v_dep = split1[1].strip()

	split2 = split1[0].split('-')
	# some of the indices won't be in the msai2seqi since the column is significant but the position on target pdb msa seq are gaps
	return '%s %s' % ('-'.join([pdbseqDict[msai2seqi[int(msai)] if int(msai) in msai2seqi else -1] for msai in split2]), v_dep)
'''

# convert
# 1042-2032-3128 0.006242240179705
# 1931-2177-3128 0.001309941125401
# 2136-3128-3140 0.003996312858620
# to
#
#
#
# protein.seqDict{} : [132 : 'B529(V)']
# msai -> seqi -> 'B529(V)'
# discard for implant alignment
'''
def sdii2resi():
	if len(sys.argv) < 5:
		print 'sdii2resi: convert msa position to residue number in pdb for a sdii result file' 
		print 'example: python utils_msa.py sdii2resi 1k2p_PF07714_full.fa.3128_3_sdii 1k2p_PF07714_full.fa 1k2p.pdb\n'
		return

	sdiifile = sys.argv[2]
	msafile = sys.argv[3]
	pdbfile = sys.argv[4]

	p = protein(pdbfile)
	m = msa(msafile)
	seqi2msai, msai2seqi = m.getPosMap(p)

	with open(sdiifile) as f:
		sdiilines = f.readlines()

	for line in sdiilines:
		print sdiiparse(line, msai2seqi, p.seqDict)
'''


def getSeqbyName():
	if len(sys.argv) < 4:
		print 'getSeqbyName: get msa sequence without gaps by searching fasta name'
		print 'example: python utils_msa.py getseqbyname PF07714_full.fa BTK_HUMAN\n'
		return

	msafile = sys.argv[2]
	msaheader = sys.argv[3].upper()
	print 'msa file: %s' % msafile
	print 'target entry: %s' % msaheader

	msaseq = ''
	m = msa(msafile)
	m.setTarget(msaheader)

	for s in m.msaArray:
		if msaheader in s[0]:
			msaheader = s[0]
			msaseq = s[1]

	outputSeq = []
	for a in msaseq:
		if a in ['.', '-', '_']:
			continue
		else:
			outputSeq.append(a)

	print msaheader
	print ''.join(outputSeq)


def getMsabyName():
	if len(sys.argv) < 4:
		print 'getMsabyName: get msa sequence with gaps by searching fasta name'
		print 'example: python utils_msa.py getmsabyname PF07714_full.fa BTK_HUMAN\n'
		return

	msafile = sys.argv[2]
	msaheader = sys.argv[3].upper()
	print 'msa file: %s' % msafile
	print 'target entry: %s' % msaheader

	msaseq = ''
	m = msa(msafile)
	m.setTarget(msaheader)

	for s in m.msaArray:
		if msaheader in s[0]:
			print s[0]
			print s[1]

def pdistDistribution():
	if len(sys.argv) < 3:
		print 'pdist: write pairwise sequence simiarity value in a file' 
		print 'example: python utils_msa.py pdist PF07714_full.fa BTK_HUMAN\n'
		return

	msafile = sys.argv[2]
	target = sys.argv[3]
	outfile = '%s.pdist' % msafile

	print 'msa file: %s\ntarget header: %s\noutfile: %s\n' % (msafile, target, outfile)

	print 'loading msa ...'
	m = msa(msafile)
	m.setTarget(target)
	print

	print 'saving to [%s] ...' % outfile
	pdist = m.getpdist()
	np.savetxt(outfile, pdist, fmt='%.8', delimiter='\n')
	print '%d pdist saved.' % len(pdist)


def reduceByHamming():
	if len(sys.argv) < 3:
		print 'reduceByHamming: reduce a msa file by selecting sequences that have sequential similarity < 62% (hamming dist > 0.38)'
		print 'example: python utils_msa.py reducebyhamming PF07714_full.fa BTK_HUMAN\n'
		return

	msafile = sys.argv[2]
	target = sys.argv[3]
	outfile = '%s.s62' % msafile

	print 'msa file: %s\ntarget header: %s\noutfile: %s\n' % (msafile, target, outfile)

	print 'loading msa ...'
	m = msa(msafile)
	m.setTarget(target)
	print 

	m.hammingReduction(outfile, 0.62)



def reduceByWeight():
	if len(sys.argv) < 5:
		print 'reduceByWeight: reduce a msa file by weighing and reduce scale (x%)'
		print 'example: python utils_msa.py reducebyweight 1k2p_PF07714_full.fa test.weight BTK_HUMAN 0.5\n'
		return

	msafile = sys.argv[2]
	weightfile = sys.argv[3]
	target = sys.argv[4]
	scale = float(sys.argv[5])
	outfile ='%s.r%d' % (msafile, scale*100)
	print 'msa file: %s' % msafile
	print 'weight file: %s' % weightfile
	print 'target: %s' % target
	print 'reduce scale: %f' % scale
	print 'output file: %s' % outfile

	weight = np.loadtxt(weightfile, delimiter=',')
	print 'weight loaded : %s' % repr(weight.shape)

	print 'loading msa file ...'
	m = msa(msafile)
	m.setTarget(target)

	rlist=[]
	for i in xrange(0, len(weight)):
		rlist.append((i, weight[i]))

	# 0 -> len(weight)
	# small -> large
	sort_rlist = sorted(rlist, key=lambda x: x[1])

	#for k in xrange(0, len(sort_rlist)):
	#	print '[%d]:[%s]' % (k, repr(sort_rlist[k]))

	goal = int(len(weight) * (1-scale))

	target_flag = False
	fout = open(outfile, 'w')
	# save msa sequences with large weights
	print 'Writing output ...'
	for k in xrange(goal, len(weight)):
		(index, w) = sort_rlist[k]
		#print '%d, %f' % (index, w)
		if m.msaArray[index][0] == m.target[0]:
			target_flag = True
		fout.write('>%s\n%s\n' % (m.msaArray[index][0], m.msaArray[index][1]))
	if target_flag == False:
		print 'Inserting target sequence: %s' % m.target[0]
		fout.write('>%s\n%s\n' % (m.target[0], m.target[1]))
	fout.close()
	print 'reduced msa: [%s]\nlen: %d' % (outfile, goal)


def MSAReduction():
	if len(sys.argv) < 4:
		print 'msareduction: reduce columns and rows by cutoffs'
		print 'example: python utils_msa.py msareduction PF07714_full.fa BTK_HUMAN 0.2 0.62\n'
		return

	msafile = sys.argv[2]
	target = sys.argv[3]
	gap_cutoff = float(sys.argv[4]) # gap cutoff
	hamming_cutoff = float(sys.argv[5]) # hamming cutoff

	print 'msa file: %s' % msafile
	print 'target: %s' % target
	print 'gap cutoff: %f' % gap_cutoff
	print 'hamming cutoff: %s' % hamming_cutoff

	print 'loading msa file ...'
	m = msa(msafile)
	m.setTarget(target)

	(score, column_index, row_index) = m.get_msaboard_RC_RR(gap_cutoff, hamming_cutoff)
	'''
	print 'score matrix:'
	for i in xrange(0, len(score)):
		print repr(score[i])
	print 'column index: %s' % repr(column_index)
	print 'row index: %s' % repr(row_index)
	'''

	#np.savetxt(msafile+'.score', score, fmt='%.8', delimiter=',')
	print 'save score to file: [%s]' % (msafile+'.score')
	np.savetxt(msafile+'.score', score, fmt='%d', delimiter=',')
	print 'save reduced indices to file: [%s]' % (msafile+'.index')
	fout = open(msafile+'.index', 'w')
	fout.write(','.join([str(i) for i in column_index])+'\n')
	fout.write(','.join([str(i) for i in row_index])+'\n')
	fout.close()






def main():

	dispatch = {
		'resi2msai': resi2msai, 'msai2resi':msai2resi, 'sdii2resi': sdii2resi, 'getseqbyname': getSeqbyName, 'getmsabyname': getMsabyName,
		'reducebyweight': reduceByWeight, 'reducebyhamming': reduceByHamming, 'resi2target': resi2target, 'pdist': pdistDistribution, 'msareduction':MSAReduction
	}

	if len(sys.argv)<2:
		for k in dispatch:
			dispatch[k]()
		return

	cmd = sys.argv[1]

	flag = False
	for key in dispatch:
		if key == cmd:
			dispatch[key]()
			flag = True
	if flag == False:
		print 'Wrong cmd string'



if __name__ == '__main__':
	main()
