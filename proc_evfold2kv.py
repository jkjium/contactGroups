import sys
from msa import msa

# convert evfold output to kv format
def main():
	if len(sys.argv) < 4:
		print 'Usage: python proc_evfold2sdii.py PF00589.dca PF00021_full.txt UPAR_MOUSE'
		print 'output: PF00589.dca.kv, PF00589.wmi.kv'
		exit()

	dcafile = sys.argv[1]
	msafile = sys.argv[2]
	seqname = sys.argv[3]

	strArray = dcafile.split('.')
	wmioutfile = strArray[0]+'.wmi.kv'
	dcaoutfile = strArray[0]+'.dca.kv'

	start_idx = '' 
	m = msa(msafile)
	posmap, msaseq = m.getPosMapbyName(seqname)
	#print '%s\n\n%s' % (repr(posmap), msaseq)

	'''
	>UPAR_MOUSE/215-294
	$ cat PF00021-UPAR_MOUSE.msa
	.............................CYSC..EGN...NTLG......CSSEea...........slINCRG.P..MNQ...................CLV.....A....
	T...G.....L...D............V.L.....G..N............RSY.TV...RGCA..T.A..SWCq........gSHVA.D.S..F.P...T.....H...L.N.
	.V.SV...................SC..CH.GSGCN.......................................
	'''	
	wmi = {}
	dca = {}
	'''
	evfold output:
	215 C 216 Y 0.350734 0.233306
	'''
	fwmi = open(wmioutfile, 'w')
	fdca = open(dcaoutfile, 'w')
	with open(dcafile) as fp:
		for line in fp:
			line = line.strip()
			evfoldArray = line.split(' ')
			key = '%d-%d' % (posmap[int(evfoldArray[0])], posmap[int(evfoldArray[2])])

			if evfoldArray[1]!=msaseq[posmap[int(evfoldArray[0])]] or evfoldArray[3]!=msaseq[posmap[int(evfoldArray[2])]]:
				print 'error: %s(%s)-%s(%s) -> %s [%s-%s]' % (evfoldArray[0], evfoldArray[1], evfoldArray[2], evfoldArray[3], key, msaseq[posmap[int(evfoldArray[0])]], msaseq[posmap[int(evfoldArray[2])]])
				exit()

			fwmi.write('%s %s\n' % (key, evfoldArray[4]))
			fdca.write('%s %s\n' % (key, evfoldArray[5]))

	print '%s.-.kv saved.' % strArray[0]
	fwmi.close()
	fdca.close()

if __name__ == '__main__':
	main()