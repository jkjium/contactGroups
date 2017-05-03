"""
Bar chart demo with pairs of bars grouped for easy comparison.
"""

'''
bar_orange = [236/256, 189/256, 137/256];
bar_red = [206/256, 147/256, 129/256];
bar_green = [191/256, 202/256, 188/256];
bar_gray = [171/256, 167/256, 168/256];
bar_brown = [238/256, 226/256, 212/256];
bar_blue = [182/256,226/256,240/256;];
'''


import numpy as np
import matplotlib.pyplot as plt


n_groups = 20

#infileA = 'cath174.b62.rmsd.stat'
#infileB = 'cath174.cb2.rmsd.stat'

infile = 'su.b62.sum.stat'

data = np.loadtxt(infile,delimiter=',')

print(repr(data))


fig, ax = plt.subplots()

index = np.arange(n_groups)
bar_width = 0.3

opacity = 0.4
error_config = {'ecolor': '0.1'}

rects1 = plt.bar(index, data[0,:], bar_width,
                 alpha=opacity,
                 color='#5293a8',
                 #yerr=stdA,
                 #error_kw=error_config,
                 label='SCSC')

rects2 = plt.bar(index + bar_width, data[1,:], bar_width,
                 alpha=opacity,
                 color='#712a14',
                 #yerr=stdB,
                 #error_kw=error_config,
                 label='BLOSUM62')

plt.xlabel('Amino acids')
plt.ylabel('Marginal substitution scores')
plt.title('Marginal ubstritution scores comparison')

penalty = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
plt.xticks(index + bar_width / 2, penalty)

plt.legend()

plt.tight_layout()
plt.show()

