import sys
import commp as cp

# generate bash script for sm generation
def cflat2sm(arglist):
	if len(arglist) < 2:
		cp._err('Usage: python utils_bash.py cflat2sm stub namestr')
	stubfile = arglist[0]	
	namestr = arglist[1]

	sarr = namestr.split(',')
	prefix = sarr[0]
	# scol
	# awk '{printf "python utils_flat.py scolsingle %s-std-mipdca.cflat tip 5 1 14 %s-tip5-d14.scol\n", $4,$4}' 3-p90-nseq1k-p90gap-pdbok-mapok.tsv
	scolcmd = []
	if len(sarr) == 4: # single: tip5-d14 5 1 14
		c = 'awk \'{printf "python utils_flat.py scolsingle %%s-std-mipdca.cflat tip %s %s %s %%s-%s.scol\\n", $1,$1}\' %s' % (sarr[1], sarr[2], sarr[3], prefix, stubfile)
	elif len(sarr) == 6: # double: tip4-m05d11 4 0 0.5 1 11
		c = 'awk \'{printf "python utils_flat.py scolinter %%s-std-mipdca.cflat tip %s %s %s %s %s %%s-%s.scol\\n", $1,$1}\' %s' % (sarr[1], sarr[2], sarr[3], sarr[4], sarr[5], prefix, stubfile)
	elif len(sarr) == 7:
		c = 'awk \'{printf "python utils_flat.py scolunion %%s-std-mipdca.cflat tip %s %s %s %s %s %%s-%s.scol\\n", $1,$1}\' %s' % (sarr[1], sarr[2], sarr[3], sarr[4], sarr[5], prefix, stubfile)

	else:
		cp._err('invalid namestr: %s' % namestr)
	scolcmd.append(c)

	# write scol sh file
	outfile = 'batch_scol.%s.psh' % prefix
	with open(outfile, 'w') as fp:
		fp.write('%s\n' % '\n'.join(scolcmd))
	cp._info('save to %s' % outfile)

	# wfreq
	# awk '{printf "python utils_pfammsa.py wfreq %s_p90.txt.score.1.flu %s_p90.txt.50.weight %s-tip5-d14.scol %s_p90.tip5-d14.50.wfreq\n", $4,$4,$4, $4}' 3-p90-nseq1k-p90gap-pdbok-mapok.tsv
	# 50, 62, 70 weight
	wfreqcmd = []
	wlist = [50, 62, 70]
	for w in wlist:
		c = 'awk \'{printf "python utils_pfammsa.py wfreq %%s_p90.txt.score.1.flu %%s_p90.txt.%d.weight %%s-%s.scol %%s_p90.%s.%d.wfreq\\n", $1,$1,$1, $1}\' %s' % (w, prefix, prefix, w, stubfile)
		wfreqcmd.append(c)
	outfile = 'batch_wfreq.%s.psh' % prefix
	with open(outfile, 'w') as fp:
		fp.write('%s\n' % '\n'.join(wfreqcmd))
	cp._info('save to %s' % outfile)

	# .allwfreq
	# cat *tip5-d14.50.wfreq > tip5-d14.50.allwfreq
	smcmd = []
	outfile = 'batch_sm.%s.psh' % prefix
	for w in wlist:
		c = 'cat *%s.%d.wfreq > %s.%d.allwfreq;python utils_pfammsa.py wfreq2sm %s.%d.allwfreq wf scsc.%s.%d' % (prefix, w, prefix, w, prefix, w, prefix, w)
		smcmd.append(c)
	with open(outfile, 'w') as fp:
		fp.write('%s\n' % '\n'.join(smcmd))
	cp._info('save %s' % outfile)
	print '\ncat batch_scol*.psh|sh > batch_scol.sh; cat batch_wfreq*.psh|sh > batch_wfreq.sh; cat batch_sm*.psh > batch_sm.sh'


def main():
	if len(sys.argv)<2:
		cp._err('Usage: python utils_bash.py cmd [args ...]')

	dispatch = {
		'cflat2sm':cflat2sm,
	}

	if sys.argv[1] in dispatch:
		dispatch[sys.argv[1]](sys.argv[2:])
	else:
		cp._err('invalid cmd: %s' % sys.argv[1])

if __name__ == '__main__':
	main()