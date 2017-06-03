import numpy as np
import matplotlib.pyplot as plt
import operator

'''
aas01 = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 
		'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

aa210 = []
for i in range(0, len(aas01)):
	for j in range(i,len(aas01)):
		aa210.append('%s-%s' % (aas01[i], aas01[j]))

print(repr(aa210))
print(len(aa210))
'''

bar_width = 0.3
opacity = 0.8
spnum = 2

h = []
#xt =[]
with open('pairfreq.list.allpairfreq.dist210.n1') as fp:
	for line in fp:
		strarr = line.strip().split(' ')
		if len(strarr)!=404:
			print('invalid line length')
		h.append(('%s-%s' % (strarr[1][0], strarr[1][1]), float(strarr[3])))

nh = np.array([v[1] for v in h])
m = np.min(nh[np.nonzero(nh)])
#print(m)

#h.sort(key=operator.itemgetter(1), reverse=True) # sort by value

data = [v-m for k,v in h]
xt = [k for k,v in h]

print(len(xt))
bar_width = 0.3
opacity = 0.8
spnum = 2

# start plot
fig ,ax = plt.subplots(spnum, figsize=(18,6), sharey=True)

#n_groups = len(freq)
index  = np.arange(105)
ax[0].bar(index, data[0:105], bar_width,
			alpha=opacity,
			color = '#75a8b9',
			label='Entropy')
plt.sca(ax[0])
plt.xticks(index+bar_width/2, xt[0:105], rotation='vertical', fontname = "monospace")
#plt.ylabel('Number of counts')
plt.legend()


index  = np.arange(105)

ax[1].bar(index, data[105:], bar_width,
			alpha=opacity,
			color = '#75a8b9',
			label='Entropy')
plt.sca(ax[1])
plt.xticks(index+bar_width/2, xt[105:], rotation='vertical', fontname = "monospace")
#plt.ylabel('Number of counts')
plt.legend()


#plt.xlabel('Pair substitution frequency')
plt.tight_layout()
plt.show()
