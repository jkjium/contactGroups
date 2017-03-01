import sys
from alignflat import alignflat

def printtotal():
	if len(sys.argv) < 3:
		print 'python utils_alignflat.py printtotal test.allflat'
		return	
	af = alignflat(sys.argv[2])
	print '%d %d %d %f %d' % (af.totalnid, af.totalnsm, af.totalngp, 1.0*af.totalnid/af.totalngp, af.totalres)


def dump():
	if len(sys.argv) < 3:
		print 'python utils_alignflat.py dump test.(all)flat'
		return
	af = alignflat(sys.argv[2])
	af.dump()


def main():

	dispatch = {
		'printtotal':printtotal, 'dump':dump
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