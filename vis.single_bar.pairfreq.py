import numpy as np
import matplotlib.pyplot as plt

n_groups = 50

xt = []
freq = []
with open('pairfreq.1.top') as fp:
	for line in fp:
		strarr = line.strip().split(' ')
		xt.append(strarr[0])
		freq.append(int(strarr[1]))

nfreq = np.array(freq)

index  = np.arange(n_groups)

bar_width = 0.3
opacity = 0.8
spnum = 2

fig ,ax = plt.subplots(spnum, figsize=(12,6), sharey=True)

ax[0].bar(index, nfreq[0:50], bar_width,
			alpha=opacity,
			color = '#75a8b9',
			label='Single change')
plt.sca(ax[0])
plt.xticks(index+bar_width/2, xt[0:50], rotation='vertical')
plt.ylabel('Number of counts')
plt.legend()

ax[1].bar(index, nfreq[50:100], bar_width,
			alpha=opacity,
			color = '#75a8b9',
			label='Single change')
plt.sca(ax[1])
plt.xticks(index+bar_width/2, xt[50:100], rotation='vertical')
plt.ylabel('Number of counts')
plt.legend()

plt.xlabel('Pair substitution frequency')
plt.tight_layout()
plt.show()


