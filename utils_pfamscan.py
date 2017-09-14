import sys
import json

'''
advanced version of proc_pfamscan_json.py

'''

class pfamscan(object):
	'''
	[
		{
		"model_length":"42",
		"align":[
				"#HMM       DvdECasgthnCpentvCvNteGsfeCvCeegye",
				"#MATCH     DvdEC  g++ C+++++C+Nt GsfeC+C +gy+",
				"#PP        9********************************8",
				"#SEQ       DVDECSLGANPCEHAGKCINTLGSFECQCLQGYT"
				],
		"env":{"to":"40","from":"2"},
		"name":"EGF_CA",
		"acc":"PF07645.13",
		"sig":1,
		"evalue":"1.2e-12",
		"desc":"Calcium-binding EGF domain",
		"hmm":{"to":"34","from":"1"},
		"act_site":null,
		"type":"Domain",
		"bits":"47.6",
		"clan":"CL0001",
		"seq":{"to":"35","from":"2","name":"A"}
		}
	]
	'''
	#def __init__(self, jsonstring):
#		jdata = json.loads(jsonstring)[0]
	def __init__(self, jdata):
		self.model_length = int(jdata['model_length'])
		aln = jdata['align']
		self.alnhmm   = str(aln[0][len('#HMM       '):])
		self.alnmatch = str(aln[1][len('#MATCH     '):])
		self.alnpp    = str(aln[2][len('#PP        '):])
		self.alnseq   = str(aln[3][len('#SEQ       '):])
		self.envfrom = int(jdata['env']['from']) - 1 
		self.envto   = int(jdata['env']['to']) - 1
		self.acc = jdata['acc']
		self.sig = int(jdata['sig'])
		self.evalue = float(jdata['evalue'])
		self.desc = jdata['desc']
		self.hmmfrom = int(jdata['hmm']['from']) - 1 
		self.hmmto   = int(jdata['hmm']['to']) - 1
		self.act_site = str(jdata['act_site'])
		self.type = jdata['type']
		self.bits = float(jdata['bits'])
		self.clan = jdata['clan']
		self.seqfrom = int(jdata['seq']['from']) - 1
		self.seqto   = int(jdata['seq']['to']) - 1
		self.seqname    = jdata['seq']['name']

		self.pfamid = self.acc[:7]



	def dump(self):
		print '\n-----------------------------'
		print 'model_length:   %s' % self.model_length
		print 'alnhmm:         %s' % self.alnhmm
		print 'alnmatch:       %s' % self.alnmatch
		print 'alnpp:          %s' % self.alnpp
		print 'alnseq:         %s' % self.alnseq
		print 'envfrom:        %d' % self.envfrom
		print 'envto:          %d' % self.envto
		print 'acc:            %s' % self.acc
		print 'pfamid:         %s' % self.pfamid
		print 'sig:            %d' % self.sig
		print 'evalue:         %e' % self.evalue
		print 'desc:           %s' % self.desc
		print 'hmmfrom:        %d' % self.hmmfrom
		print 'hmmto:          %d' % self.hmmto
		print 'act_site:       %s' % self.act_site
		print 'type:           %s' % self.type
		print 'bits:           %f' % self.bits
		print 'clan:           %s' % self.clan
		print 'seqfrom:        %d' % self.seqfrom
		print 'seqto:          %d' % self.seqto
		print 'seqname:        %s' % self.seqname
		print '-----------------------------\n'


	# write matched HMM to as fasta file
	def writeHMMfa(self, outfile):
		with open(outfile, 'w') as fp:
			fp.write('>%s\n%s' % (self.seqname ,self.alnhmm.upper()))


# one json file may contain multiple matches
# parse all the matches into a list
class utils_pfamscan(object):
	#use generator!!
	def __init__(self, jsonfile):
		self.name = jsonfile
		with open(jsonfile) as fp:
			self.pslist = [pfamscan(js) for js in json.loads(fp.read())] # empty line does not matter

	def dump(self):
		for ps in self.pslist:
			ps.dump()
		print '%d pfamscan records printed' % len(self.pslist)

	# return the first PFXXXXX match
	def getMatchpfs(self, pfamid):
		for ps in self.pslist:
			if ps.pfamid[0:7] == pfamid:
				return ps
		return False

	# get hmm sequence in json
	def getHmmSeq(self, opt='first'):
		if opt == 'first':
			hmmfa = self.pslist[0].alnhmm
		elif 'PF' in opt and len(opt) == 7:
			ps = self.getMatchfs(opt)
			if ps == False:
				cp._err('%s not found.' % opt)
			hmmfa = ps.alnhmm
		elif opt == 'all':
			hmmfa = '\n'.join([ '>%s_hmm\n%s' % (ps.seqname, ps.alnhmm) for ps in self.pslist])


def test():

	ups = utils_pfamscan('t.json')
	print 'json loaded: -----------------------'
	ups.dump()
	print 'get PF02576 match -----------------------'
	ps = ups.getMatchpfs('PF02576')
	ps.dump()



	# test writeHMMfa
	'''
	with open('t.json') as fp:
		jsonstr = fp.readline()
	pfs = pfamscan(jsonstr)
	pfs.dump()
	'''
	#outfile = 't.hmmfa'
	#pfs.writeHMMfa(outfile)
	#print 'file %s saved.' % outfile



# main routine
def main():
	if len(sys.argv)<2:
		print 'Usage: python utils_protein.py cmd pdbfile [args ...]'
		return

	dispatch = {
		'test':test
	}

	cmd = sys.argv[1]

	if sys.argv[1] not in dispatch:
		print 'invalid cmd: %s' % sys.argv[1]
		return
	else:
		dispatch[sys.argv[1]]()





if __name__ == '__main__':
	main()