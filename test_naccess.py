'''
test naccess class 
'''

from naccess import naccess
from naccess import rsa

def main():
	rsafile = '1k2p.rsa'
	na = naccess(rsafile)
	na.dump()

if __name__ == '__main__':
	main()
