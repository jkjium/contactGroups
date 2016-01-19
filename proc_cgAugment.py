import sys
from cgroup import cgroup

def main():
	if len(sys.argv) < 3:
		print 'Procedure for adding additional properties on contact groups'
		print 'Usage: python proc_cgAugment.py AAtips.def cgfile'
		print '\tAAtips.def: Amino Acid property description file'
		print '\tcgfile: contact group file, each line is a contact group'
		print '\tcontact group example: 2nrt.tip,LEY,380 364 417'
		return

	AADesFile = sys.argv[1] 
	print AADesFile
	CGFile = sys.argv[2]
	print CGFile

	# split in lines without '\n'
	print 'loading AA description ...'
	with open(AADesFile) as f:
		AAlines = f.read().splitlines()

	# compose type map from Amino Acid property file
	# AAtips.def
	typeMap={}
	for line in AAlines:
		strArr = line.split(',')
		typeMap[strArr[1]]=strArr[3]

	fout = open(CGFile+'.aug', 'w')
	# load contact groups
	with open(CGFile) as f:
		CGlines = f.read().splitlines()
	# output augmented contact group strings
	print 'adding properties ...'
	for line in CGlines:
		cg = cgroup(line)
		cg.getType(typeMap)
		fout.write(cg.getString())

	fout.close()
	print 'done. output: ' + CGFile+'.aug'

	'''
	cgStr = '2nrt.tip,LEY,380 364 417'
	cg = cgroup(cgStr)
	print cg.getString()

	cg.getType(typeMap)
	print cg.getString()
	'''


if __name__=='__main__':
	main()