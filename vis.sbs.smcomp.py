import numpy as np
import matplotlib.pyplot as plt

b62 = np.loadtxt('B62')
scsc = np.loadtxt('scsc')
scsc23 = np.loadtxt('scsc2.3')
outname = 'sm_comp.png'

AA = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']

xt =[]
b62a = []
scsca = []
scsc23a = []
for i in range(0, len(AA)):
	for j in range(i, len(AA)):
		xt.append('%s-%s' % (AA[i], AA[j]))
		b62a.append(b62[i][j])
		scsca.append(scsc[i][j])
		scsc23a.append(scsc23[i][j])

print(len(xt))
n_groups=itv=55
spnum = 4
fig, ax = plt.subplots(spnum, figsize=(18,11), sharey=True)

bar_width = 0.25
opacity = 1.0

for i in range(0, spnum):
		index = np.arange(n_groups if (i+1)*itv <210 else (210 - i*itv))
		sxt = xt[i*itv:((i+1)*itv if (i+1)*itv <210 else 210)]

		ax[i].bar(index, b62a[i*itv:((i+1)*itv if (i+1)*itv <210 else 210)], bar_width,
						alpha=opacity,
						color='#75a8b9',
						label='B62')

		ax[i].bar(index + bar_width, scsca[i*itv:((i+1)*itv if (i+1)*itv <210 else 210)], bar_width,
						alpha=opacity,
						color='#8d5543',
						label='SCSC')

		ax[i].bar(index + 2*bar_width, scsc23a[i*itv:((i+1)*itv if (i+1)*itv <210 else 210)], bar_width,
						alpha=opacity,
						color='#f3db81',
						label='SCSC2')
		plt.sca(ax[i])
		plt.xticks(index + bar_width / 2, sxt, rotation='vertical')

		plt.ylabel('Log-odds scores')
		plt.legend()

plt.xlabel('Amino acid pairs',fontsize='14')        

plt.tight_layout()
#plt.suptitle('%s residue contacts' % group.title(),fontsize='14')
plt.show()

fig.savefig(outname)
