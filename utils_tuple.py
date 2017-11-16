import sys
import numpy as np
import commp as cp

from utils_mprun import qrun

# data: np matrix
# varset: variable tuple
def hat(data, varset): 
	nhatlist = []
	nhatmat = []
	for t in data[:,varset]:
		aatuple =  [cp.scoreaa['aa'][int(i)] for i in t]
		nhat = [cp.aaprop[a][21] for a in aatuple]
		nhatlist.append(sum(nhat))
		nhatmat.append(nhat)
	# list of sum hat per tuple
	npnhat = np.array(nhatlist)
	# nx3 (order) matrix
	npmhat = np.array(nhatmat)
	#cp._info(repr(nhatmat))
	eig = np.linalg.eig(np.dot(npmhat.T,npmhat))
	return '%s,%.8f,%.8f\n' % (','.join([str(int(i)) for i in varset]), npnhat.std()/npnhat.mean(), max(eig[0]))

# single argument to fit qrun scheme
def mphat(args):
	return hatcv(args[0], args[1])

# calculate Coefficient of variation for each tuple column
def tuplehat(arglist):
	if len(arglist) < 5:
		cp._err('Usage: python utils_tuple.py mp tuplehat scorefile rcolfile order outfile')

	opt = arglist[0]
	if opt not in ['mp', 'sp']:
		cp._err('invalid opt %s' % opt)

	scorefile = arglist[1]	# score file by category
	rcolfile = arglist[2]	# columns involved in the significance calculation
	order = int(arglist[3])
	outfile = arglist[4]

	score = np.loadtxt(scorefile, delimiter=',')
	# one line: 1,3,5,10 ...
	if rcolfile == 'na':
		rcol = np.array([i for i in xrange(0, score.shape[1])])
	else:
		rcol = np.loadtxt(rcolfile, delimiter=',')
		rcol -= 1


	# generate 
	tasks = cp.ncrvar(rcol, order) 

	if opt == 'sp': # single process version
		with open(outfile, 'w') as fp:
			for t in tasks:
				fp.write(hat(score, t))
	elif opt == 'mp':
		mptasks = [(score, t) for t in tasks]
		qrun(mphat, mptasks, outfile, 3)

	cp._info('save to %s' % outfile)


#
def main():
	if len(sys.argv)<2:
		cp._err('Usage: python utils_tuple.py cmd [args ...]')

	dispatch = {
		'tuplehat':tuplehat
	}

	if sys.argv[1] in dispatch:
		dispatch[sys.argv[1]](sys.argv[2:])
	else:
		cp._err('invalid cmd: %s' % sys.argv[1])

if __name__ == '__main__':
	main()