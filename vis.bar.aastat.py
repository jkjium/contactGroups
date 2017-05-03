"""
Bar chart demo with pairs of bars grouped for easy comparison.
"""

import numpy as np
import matplotlib.pyplot as plt

'''
xt = ['G','P','C','U','A','I','L','V','M','F','W','N','Q','S','T','Y','D','E','R','H','K']
n_groups = len(xt)


infile = '13-rpdb.list.aa.stat'
outfile = infile+'.png'

data = np.loadtxt(infile,delimiter=' ')

fig, ax = plt.subplots(figsize=(8,4))

index = np.arange(n_groups)
bar_width = 0.3
opacity = 0.4
#error_config = {'ecolor': '0.1'}
rects1 = plt.bar(index, data, bar_width,
                 alpha=opacity,
                 color='#5293a8',
                 #yerr=stdA,
                 #error_kw=error_config,
                 #label='SCSC')
                )

plt.xlabel('Amino acid')
plt.ylabel('Frequency')
plt.title('Amino acid frequency')
plt.xticks(index + bar_width / 2, xt)
#plt.legend()

plt.tight_layout()
plt.show()

fig.savefig(outfile)
'''
