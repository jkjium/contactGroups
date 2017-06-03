import numpy as np
import matplotlib.pyplot as plt


aas01 = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 
		'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
xt =[]

for a in aas01:
	for b in aas01:
		xt.append('%s-%s' % (a,b))
#print(repr(xt))
print(len(xt))

bar_width = 0.3
opacity = 0.8
spnum = 4


outtitle = 'n3'
infile = 'pairfreq.list.allpairfreq.dist210.'+outtitle

with open(infile) as fp:
	for line in fp:
		strarr = line.strip().split(' ')
		if len(strarr)!=404:
			print('invalid line length')
			continue
		outfile = 'pf.%03d_%s_%s.png' % (int(strarr[0]), strarr[1],outtitle )
		print(outfile)
		nfreq = np.array([float(v) for v in strarr[4:]])

		n_groups = 100
		index  = np.arange(n_groups)

		fig ,ax = plt.subplots(spnum, figsize=(18,9), sharey=True)

		ax[0].bar(index, nfreq[0:100], bar_width,
					alpha=opacity,
					color = '#75a8b9',
					label='Single change')
		plt.sca(ax[0])
		plt.xticks(index+bar_width/2, xt[0:100], rotation='vertical', fontname = "monospace")
		plt.ylabel('Normalized frequency')
		#plt.legend()

		ax[1].bar(index, nfreq[100:200], bar_width,
					alpha=opacity,
					color = '#75a8b9',
					label='Single change')
		plt.sca(ax[1])
		plt.xticks(index+bar_width/2, xt[100:200], rotation='vertical', fontname = "monospace")
		plt.ylabel('Normalized frequency')
		#plt.legend()
		ax[2].bar(index, nfreq[200:300], bar_width,
					alpha=opacity,
					color = '#75a8b9',
					label='Single change')
		plt.sca(ax[2])
		plt.xticks(index+bar_width/2, xt[200:300], rotation='vertical', fontname = "monospace")
		plt.ylabel('Normalized frequency')
		#plt.legend()


		ax[3].bar(index, nfreq[300:], bar_width,
					alpha=opacity,
					color = '#75a8b9',
					label='Single change')
		plt.sca(ax[3])
		plt.xticks(index+bar_width/2, xt[300:], rotation='vertical', fontname = "monospace")
		plt.ylabel('Normalized frequency')
		#plt.legend()		

		plt.xlabel('[%s-%s]: Normalized freq: %.4f, Entropy: %.4f' % (strarr[1][0], strarr[1][1], float(strarr[2]), float(strarr[3])))
		plt.tight_layout()
		#plt.show()
		fig.savefig(outfile)
		fig.clf()
		plt.close()
		
		#break
		


