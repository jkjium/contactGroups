import numpy as np
import matplotlib.pyplot as plt

N = 5
b62 = (18.61,8.91)
#scsc = (2, 3, 4, 1, 2)
scsc = (27.19,17.11)
N = len(b62)
#omen_std = (3, 5, 2, 3, 3)

ind = np.arange(N)  # the x locations for the groups
width = 0.35       # the width of the bars
opacity = 1.0

fig, ax = plt.subplots()
rects1 = ax.bar(ind, b62, width, alpha=opacity, color = '#75a8b9')
rects2 = ax.bar(ind + width, scsc, width, alpha=opacity, color='#f3db81')

# add some text for labels, title and axes ticks
ax.set_ylabel('Ratio')
ax.set_title('True / False positive ratio by Matrix')
ax.set_xticks(ind + width / 2)
ax.set_xticklabels(('Total sum rate', 'Case-wised rate'))

ax.legend((rects1[0], rects2[0]), ('B62', 'SCSC'))


def autolabel(rects):
    """
    Attach a text label above each bar displaying its height
    """
    for rect in rects:
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width()/2., 0.9*height,
                '%.2f' % float(height),
                ha='center', va='bottom')

autolabel(rects1)
autolabel(rects2)

plt.show()