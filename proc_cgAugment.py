import sys
from cgroup import cgroup

def main():
	cgStr = '2nrt.tip,LEY,380 364 417'
	cg = cgroup(cgStr)
	print cg.getString()

	# split in lines without '\n'
	with open('AAtips.def') as f:
		AAlines = f.read().splitlines()
	print AAlines

	typeMap={}
	for line in AAlines:
		strArr = line.split(',')
		typeMap[strArr[1]]=strArr[3]
	print typeMap

	cg.getType(typeMap)
	print cg.getString()


if __name__=='__main__':
	main()