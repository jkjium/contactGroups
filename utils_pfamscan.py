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
	def __init__(self, jdata):
		self.model_length = int(jdata['model_length'])
		aln = jdata['align']
		self.alnhmm   = aln[0][len('#HMM       '):]
		self.alnmatch = aln[1][len('#MATCH     '):]
		self.alnpp    = aln[2][len('#PP        '):]
		self.alnseq   = aln[3][len('#SEQ       '):]
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


class utils_pfamscan(object):
	def __init__(self, jsonfile):
		with open(jsonfile) as fp:
			line = fp.readline()
		self.pslist = [pfamscan(js) for js in json.loads(line.strip())]

	def dump(self):
		for ps in self.pslist:
			ps.dump()
		print '%d pfamscan records printed' % len(self.pslist)


def main():
	ups = utils_pfamscan('t.json')
	ups.dump()

if __name__ == '__main__':
	main()