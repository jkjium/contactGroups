'''
get resid list of varname
'''
import sys
from naccess import naccess
from naccess import rsa

def main():
	if len(sys.argv)<3:
		print 'Usage: utils_naccess.py 1k2p.rsa DB'
		return

	rsafile = sys.argv[1]
	varname = sys.argv[2]
	na = naccess(rsafile)
	#na.dump()
	#na.dumpResiMap()
	print '[%s]: %s' % (varname, na.getResiList(varname))

if __name__ == '__main__':
	main()
