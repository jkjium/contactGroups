"""
Bar chart demo with pairs of bars grouped for easy comparison.
"""

'''
rcg.nearest.stat
0 pdbs.ca.nearest.2.cg
1 pdbs.sgc.nearest.2.cg
2 pdbs.tip.nearest.2.cg


rcg.cutoff.stat
0 pdbs.ca.cutoff.6.5_6.cg
1pdbs.sgc.cutoff.6.5_6.cg
2 pdbs.tip.cutoff.6.5_6.cg

'''


import numpy as np
import matplotlib.pyplot as plt



group = 'nearest'
group = 'cutoff'

infile = 'rcg.%s.stat' % group
outname = 'rcg.%s.stat.png' % group


stat = np.loadtxt(infile,delimiter=' ')


n_groups=itv=60

AA = [
    'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 
    'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'U','V', 'W', 'Y',
     ]
# init xticks
xt = []
for i in range(0, len(AA)):
        for j in range(i, len(AA)):
                xt.append('%s-%s' % (AA[i], AA[j]))

print(len(xt))


spnum = 4

fig, ax = plt.subplots(spnum, figsize=(18,11), sharey=True)


#index = np.arange(n_groups)
bar_width = 0.25
opacity = 1.0
#error_config = {'ecolor': '0.1'}

for i in range(0, spnum):
        index = np.arange(60 if (i+1)*itv <231 else (231 - i*itv))
        sxt = xt[i*60:(i*60+60 if i*60+60 <231 else 231)]

        ax[i].bar(index, stat[0,i*itv:((i+1)*itv if (i+1)*itv <231 else 231)], bar_width,
                         alpha=opacity,
                         color='#75a8b9',
                         #yerr=stdA,
                         #error_kw=error_config,
                         label='CA')

        ax[i].bar(index + bar_width, stat[1,i*itv:((i+1)*itv if (i+1)*itv <231 else 231)], bar_width,
                         alpha=opacity,
                         color='#8d5543',
                         #yerr=stdB,
                         #error_kw=error_config,
                         label='SC Centroid')

        ax[i].bar(index + 2*bar_width, stat[2,i*itv:((i+1)*itv if (i+1)*itv <231 else 231)], bar_width,
                         alpha=opacity,
                         color='#f3db81',
                         #yerr=stdB,
                         #error_kw=error_config,
                         label='Tip Atom')
        plt.sca(ax[i])
        #plt.ylim(0,9000)
        plt.xticks(index + bar_width / 2, sxt, rotation='vertical')

        plt.ylabel('Total number of contacts')
        plt.legend()

plt.xlabel('Amino acid pairs',fontsize='14')        
#fig.text(0.04, 0.5, 'common Y', va='center', rotation='vertical')


plt.tight_layout()
plt.suptitle('%s residue contacts' % group.title(),fontsize='14')
plt.show()

fig.savefig(outname)
