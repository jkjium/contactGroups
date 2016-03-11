'''
get msa position id 
'''
import sys
from msa import msa
from protein import protein

# given a residue number output the corresponding position in msa file
def resi2msai():
	if len(sys.argv) < 5:
		print 'r2p: given a residue number output the corresponding position in msa'
		print 'resi2msai()::python utils_msa.py r2p 1k2p_PF07714_full.fa 1k2p.pdb B641'
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


def msai2resi():
	pass




# function for parsing sdii result
# msai -> seqi -> 'B529(V)'
# pdbseqDict: 132 : 'B529(V)'
# sdiline: [1042-2032-3128 0.006242240179705]
def sdiiparse(sdiiline, msai2seqi, pdbseqDict):
	split1 = sdiiline.split(' ')
	v_dep = split1[1].strip()

	split2 = split1[0].split('-')
	# some of the indices won't be in the msai2seqi since the column is significant but the position on target pdb msa seq are gaps
	return '%s %s' % ('-'.join([pdbseqDict[msai2seqi[int(msai)] if int(msai) in msai2seqi else -1] for msai in split2]), v_dep)


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
def sdii2resi():
	if len(sys.argv) < 5:
		print 's2r:	convert msa position to residue number in pdb for a sdii result file' 
		print 'sdii2resi()::python utils_msa.py s2r 1k2p_PF07714_full.fa.3128_3_sdii 1k2p_PF07714_full.fa 1k2p.pdb'
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



def getSeqbyName():
	if len(sys.argv) < 4:
		print 'getSeqbyName: get msa sequence with fasta name'
		print 'getSeqbyName():: python utils_msa.py getSeqbyName PF07714_full.fa BTK_HUMAN'
		return

	msafile = sys.argv[2]
	msaheader = sys.argv[3]
	m = msa(msafile)
	for s in m.msaArray:
		if msaheader.upper() == s[0].upper():
			print s[0]
			print s[1]



def main():

	dispatch = {
		'r2p': resi2msai, 'p2r':msai2resi,
		's2r': sdii2resi, 'getSeqbyName': getSeqbyName
	}

	if len(sys.argv)<2:
		print 'Usage: utils_msa.py cmd ...'
		return

	cmd = sys.argv[1]

	for key in dispatch:
		if key == cmd:
			dispatch[key]()

	'''
	if len(sys.argv)<4:
		print 'Usage: utils_msa.py msafile pdbfile chain+resi'
		return
	msafile = sys.argv[1]
	pdbfile = sys.argv[2]
	resi = sys.argv[3]

	m = msa(msafile)
	p = protein(pdbfile)
	#print p.seq
	#print p.resDict
	resIdx = p.resDict[resi]
	posMap = m.getPosMap(p)
	#print m.msaArray[0][1]
	for i in posMap:
		if i == resIdx[0]:
			print '[Res: %s] : (seqi: %d (%s) - msai: %d (%s))' % (resi, resIdx[0], resIdx[1], posMap[i], m.msaArray[0][1][posMap[i]])
			break
	'''

if __name__ == '__main__':
	main()
