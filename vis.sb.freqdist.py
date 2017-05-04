import numpy as np
import matplotlib.pyplot as plt


aas01 = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 
		'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

aa210 = []
for i in range(0, len(aas01)):
	for j in range(i,len(aas01)):
		aa210.append('%s-%s' % (aas01[i], aas01[j]))

print(repr(aa210))
print(len(aa210))

bar_width = 0.3
opacity = 0.8
spnum = 2

freq = []
h = []
xt =[]
with open('pairfreq.list.allpairfreq.dist210') as fp:
	for line in fp:
		strarr = line.strip().split(' ')
		if len(strarr)!=404:
			print('invalid line length')
		xt.append('%s-%s' % (strarr[1][0], strarr[1][1]))
		freq.append(float(strarr[2]))
		h.append(float(strarr[3]))

diff = list(set(aa210) - set(xt))
print(repr(diff))
nfreq = np.array(freq)
th = np.array(h)
nh = th - np.min(th[np.nonzero(th)])

bar_width = 0.3
opacity = 0.8
spnum = 2

data = nfreq

fig ,ax = plt.subplots(spnum, figsize=(18,6), sharey=True)

#n_groups = len(freq)
index  = np.arange(105)
ax[0].bar(index, data[0:105], bar_width,
			alpha=opacity,
			color = '#75a8b9',
			label='Normalized frequency')
plt.sca(ax[0])
plt.xticks(index+bar_width/2, xt[0:105], rotation='vertical', fontname = "monospace")
#plt.ylabel('Number of counts')
plt.legend()


#index  = np.arange(104)
ax[1].bar(index, data[105:], bar_width,
			alpha=opacity,
			color = '#75a8b9',
			label='Normalized frequency')
plt.sca(ax[1])
plt.xticks(index+bar_width/2, xt[105:], rotation='vertical', fontname = "monospace")
#plt.ylabel('Number of counts')
plt.legend()


#plt.xlabel('Pair substitution frequency')
plt.tight_layout()
plt.show()
