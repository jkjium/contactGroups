import commp as cp
import numpy as np
import sys


def _sbatchstr(name, cmd, hours):
	sbatch=[]
	sbatch.append('#!/bin/bash')
	sbatch.append('#SBATCH --nodes=1')
	sbatch.append('#SBATCH --cpus-per-task=8')
	sbatch.append('#SBATCH --mem=64G')
	
	sbatch.append('\n#SBATCH --time=%d:0:0\n' % hours)
	print('%s %d' % (name, hours))

	sbatch.append('#SBATCH --output=job.%s.out' % name)
	sbatch.append('#SBATCH --error=job.%s.err\n' % name)

	# module
	sbatch.append('module load py-numpy/1.15.2-py2-fdkji5s')
	sbatch.append('module load py-scipy/1.1.0-py2-rrg4qpb\n')

	# cmd
	sbatch.append(cmd)

	return sbatch

def sbatchfiles(args):
	assert len(args)==2, 'Usage: python sbatch2.py sbatchfiles job.stub hours'
	hours = int(args[1])
	for j in cp.loadtuples(args[0], delimiter=','): # name, cmd
		name = j[0]
		cmd = j[1]
		sbatchstr = _sbatchstr(name, cmd, hours)
		outfile = '%s.sh' % name
		with open(outfile, 'w') as fout:
			fout.write('%s\n\n' % ('\n'.join(sbatchstr)))
		cp._info('save to %s\n' % outfile)
	

if __name__=='__main__':
	cp.dispatch(__name__)
