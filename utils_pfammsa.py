import sys
import commp as cp

"""
updated utils_msa.py
follow cls -> utils_cls => proc_
"""
class pfammsa(object):
	"""
	data structure of pfam MSA

	"""
	def __init__(self, msafile):
		self.headdict = {}
		self.msalist = []
		count = 0
		for head, seq in cp.fasta_iter(msafile):
			#print '%d\n%s\n%s\n' % (count, head, seq)
			self.msalist.insert(count, seq)
			self.headdict[head] = count
			count+=1
		self.msalen = len(self.msalist[0])
		self.msanum = len(self.msalist)

	def dump(self):
		for k in self.headdict:
			print ">%s\n%s\n" % (k, self.msalist[self.headdict[k]])


class utils_pfammsa(object):
	"""
	general functions operated on pfammsa

	"""
	def __init__(self, msafile):
		self.pm = pfammsa(msafile)

	def dump(self):
		self.pm.dump()

	def fetchbyname(self, name):
		# >H2U8X0_TAKRU/9-86
		return ('>%s\n%s\n' % (name, self.pm.msalist[self.pm.headdict[name]])) if name in self.pm.headdict else False


def main():
	up = utils_pfammsa('PF00595_test.txt')
	up.dump()
	print up.fetchbyname('F7E5S2_MACMU/1044-1120')

if __name__ == '__main__':
	main()