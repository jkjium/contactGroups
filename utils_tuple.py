import sys
import numpy as np
import commp as cp

from utils_mprun import qrun

# calculate Coefficient of variation for each tuple column
def transscore(arglist):
	if len(arglist) < 4:
		cp._err('Usage: python utils_tuple.py transscore PF00098_p90.txt.score aa 21 outscorefile')

	scorefile = arglist[0]
	srcname = arglist[1]
	dstnameidx = int(arglist[2])
	outfile = arglist[3]

	score = np.loadtxt(scorefile, delimiter=',')
	fp = open(outfile, 'w')
	for i in xrange(0, score.shape[0]):
		aaline = [cp.scoreaa[srcname][int(k)] for k in score[i,:]]
		dstscoreline = [str(cp.aaprop[a][dstnameidx]) for a in aaline]
		fp.write('%s\n' % ','.join(dstscoreline))
	fp.close()
	cp._info('save to %s' % outfile)


# data: np matrix
# varset: variable tuple
def hat(data, varset, flag): 
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
	cp._info(repr(varset))
	cp._info(repr(nhatmat))
	eig = np.linalg.eig(np.dot(npmhat.T,npmhat))
	# variable triplet,sig,mean,std,cv,max 3x3matrix eig
	return '%s,%d,%.8f,%.8f,%.8f,%.8f\n' % (
				','.join([str(int(i)) for i in varset]), 
				flag,
				npnhat.mean(), 
				npnhat.std(), 
				npnhat.std()/npnhat.mean(), 
				max(eig[0])
			)


# single argument to fit qrun scheme
def mphat(args):
	return hat(args[0], args[1], args[2])

# calculate Coefficient of variation for each tuple column
def tupleppt(arglist):
	if len(arglist) < 6:
		cp._err('Usage: python utils_tuple.py mp tupleppt scorefile rcolfile stuplefile order outfile')

	opt = arglist[0]
	if opt not in ['mp', 'sp']:
		cp._err('invalid opt %s' % opt)

	scorefile = arglist[1]	# score file by category
	rcolfile = arglist[2]	# columns involved in the significance calculation
	stuplefile = arglist[3]
	order = int(arglist[4])
	outfile = arglist[5]

	score = np.loadtxt(scorefile, delimiter=',')
	# score column index, one line: 1,3,5,10 ...
	if rcolfile == 'na':
		rcol = np.array([i for i in xrange(0, score.shape[1])])
	else:
		rcol = np.loadtxt(rcolfile, delimiter=',')
		rcol -= 1

	# significant sequential id of column index
	stuple = np.loadtxt(stuplefile, delimiter=',')
	stuple -= 1

	# convert to score column index
	stuple_set = [set(rcol[t]) for t in stuple.tolist()]

	# generate 
	tasks = cp.ncrvar(rcol, order) 

	if opt == 'sp': # single process version
		count = 0
		with open(outfile, 'w') as fp:
			for t in tasks:
				count+=1
				if count%100==0:
					print (count, repr(t))
				fp.write(hat(score, t, int(set(t) in stuple_set)))
		cp._info('total %d tuples processed' % count)
	elif opt == 'mp':
		mptasks = [(score, t, int(set(t) in stuple_set)) for t in tasks]
		qrun(mphat, mptasks, outfile, 3)
	cp._info('save to %s' % outfile)


#
def main():
	if len(sys.argv)<2:
		cp._err('Usage: python utils_tuple.py cmd [args ...]')

	dispatch = {
		'tupleppt':tupleppt,
		'transscore': transscore
	}

	if sys.argv[1] in dispatch:
		dispatch[sys.argv[1]](sys.argv[2:])
	else:
		cp._err('invalid cmd: %s' % sys.argv[1])

if __name__ == '__main__':
	main()