import sys
import numpy as np
# informative sdii filter 
# from shadpw algorithm

def main():
	if len(sys.argv) < 4:
		print 'Usage: python proc_topsdii.py PF00001_full.txt.2_sdii int(cutoff level)=3 colidx'
		print 'output: PF00001_full.txt.2_sdii.top'
		exit()

	sdiifile = sys.argv[1]
	d = int(sys.argv[2])
	col = int(sys.argv[3])
	outfile = sdiifile+('.%d.top' % d)

	# Array of tuple (2-3, 0.25)
	sdiiArray = []
	sdiiValue = []
	with open(sdiifile) as fp:
		for line in fp:
			line = line.strip()
			if len(line) < 1:
				print 'error sdii line: %s' % line
			valueArray = line.split(' ')
			sdiiArray.append((valueArray[0], float(valueArray[col])))
			sdiiValue.append(float(valueArray[col]))

	#for a in sdiiArray:
	#	print repr(a)
	#print 'sdiiArray: %s' % (repr(sdiiArray))
	#print 'sdiiValue: %s' % (repr(sdiiValue))

	sdiinp = np.array(sdiiValue)
	outlier = sdiinp.mean()+sdiinp.std()
	#print 'm: %.4f, s: %.4f, o: %.4f' % (sdiinp.mean(), sdiinp.std(), outlier)

	sdii_no_outlier = [v for v in sdiiValue if v < outlier]
	sdiinp = np.array(sdii_no_outlier)
	cutoff = sdiinp.mean() + d*sdiinp.std()
	#print 'sdii_no_outlier: %s' % repr(sdii_no_outlier)
	#print 'm1: %.4f, s1: %.4f, cutoff: %.4f' % (sdiinp.mean(), sdiinp.std(), cutoff)

	c = 0
	fout = open(outfile, 'w')
	for (var, value) in sdiiArray:
		if value > cutoff:
			fout.write('%s %.8f\n' % (var, value))
			c+=1
	fout.close()

	print '%s: meanMI: %.4f cutoff: %.4f #ofIPV: %d/%d' % (outfile, sdiinp.mean(), cutoff, c, len(sdiiValue))

if __name__ == '__main__':
	main()