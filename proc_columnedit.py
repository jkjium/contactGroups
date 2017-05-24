import sys
from collections import defaultdict
"""
For column modifications 
combine (union), intersection, substraction

"""

def main():
	if len(sys.argv) < 4 or sys.argv[2] not in ['+','-','.']: 
		print 'Usage: python proc_columnedit.py column1.txt operator column2.txt'
		print '+: union\n-: substraction\n.:intersection'
		return

	col1file = sys.argv[1]
	op = sys.argv[2]
	col2file = sys.argv[3]

	with open(col1file) as fp:
		col1 = set(fp.readlines())

	with open(col2file) as fp:
		col2 = set(fp.readlines())

	if op == '+':
		outcol = col1.union(col2)
	elif op == '-':
		outcol = col1 - col2
	elif op == '.':
		outcol = col1.intersection(col2)

	print ''.join(sorted(list(outcol)))

if __name__ == '__main__':
	main()