import sys
import commp as cp
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram
from scipy.cluster.hierarchy import linkage

def barplot_fgdist(arglist):
    if len(arglist) < 1:
        cp._err('Usage: python utils_vis_sm.py barplot_fgdist data.tsv outfile')
    infile = arglist[0]

    bar_width = 0.3
    opacity = 0.8
    spnum = 2
    title = ''

    data = np.loadtxt(infile, delimiter=' ')
    # be consistent with utils_pfammsa.py::wfreqcsfgdist
    xt = ['%s-%s' % (cp.aas01[i], cp.aas01[j]) for i in xrange(len(cp.aas01)) for j in xrange(i,len(cp.aas01))]
    fig ,ax = plt.subplots(spnum, figsize=(18,6), sharey=True)

    index  = np.arange(105)
    ax[0].bar(index, data[0:105], bar_width,
                alpha=opacity,
                color = '#75a8b9',
                label=title)
    plt.sca(ax[0])
    plt.xticks(index+bar_width/2, xt[0:105], rotation='vertical', fontname = "monospace")
    #plt.ylabel('Number of counts')
    plt.legend()

    ax[1].bar(index, data[105:], bar_width,
                alpha=opacity,
                color = '#75a8b9',
                label=title)
    plt.sca(ax[1])
    plt.xticks(index+bar_width/2, xt[105:], rotation='vertical', fontname = "monospace")
    #plt.ylabel('Number of counts')
    plt.legend()
    plt.title(infile)

    #plt.xlabel('Pair substitution frequency')
    plt.tight_layout()
    #plt.show()
    fig.savefig(infile+'.png')
    cp._info('save to %s.png' % infile)


def barplot_bgdist(arglist):
    if len(arglist) < 2:
        cp._err('Usage: python utils_vis_sm.py barplot_bgdist data.tsv outfile')
    infile = arglist[0]
    outfile = arglist[1]

    data = np.loadtxt(infile,delimiter=' ')
    n = len(data)
    fix, ax = plt.subplots()
    index = np.arange(n)

    bar_width = 0.3
    opacity = 0.4

    rects = plt.bar(index+bar_width / 2, data, bar_width,
                    alpha=opacity,
                    color='#5293a8',
                    label='background distribution')
    #plt.xlabel('Amino acids')
    #plt.ylabel('')
    #plt.title('Marginal ubstritution scores comparison')
    plt.ylim(0,0.2)
    xlabels = cp.aat01
    plt.xticks(index + bar_width / 2, xlabels)

    plt.legend()

    plt.tight_layout()
    #plt.show()
    plt.savefig(outfile)
    cp._info('save to %s' % outfile)

def _augmented_dendrogram(*args, **kwargs):
    ddata=dendrogram(*args, **kwargs)
    if not kwargs.get('no_plot', False):
        for i,d in zip(ddata['icoord'], ddata['dcoord']):
            x=0.5*sum(i[1:3])
            y=d[1]
            plt.plot(x,y,'ro')
            str = '%.3g' % y
            plt.annotate(str, (x,y), xytext=(0,-8), textcoords='offset points', va='top', ha='center')
    return ddata

def dendrogram_fgdist(arglist):
    if len(arglist) < 1:
	    cp._err('Usage: python utils_vis_sm.py dendrogram_fgdist cs3.210dist.vec')
    sys.setrecursionlimit(10000)

    datafile = arglist[0]
    outfile = '%s.dendrogram.png' % datafile
    x = np.loadtxt(datafile, delimiter=' ')

    linkage_matrix = linkage(x, "single", metric='correlation')
    plt.figure(figsize=(24, 15))
    #ddata = dendrogram(linkage_matrix)
    xt = ['%s-%s' % (cp.aas01[i], cp.aas01[j]) for i in xrange(len(cp.aas01)) for j in xrange(i,len(cp.aas01))]
    ddata = dendrogram(linkage_matrix, labels=xt, leaf_rotation=90, color_threshold=1)
    #ddata = _augmented_dendrogram(linkage_matrix, labels=xt, leaf_rotation=90, color_threshold=1)
    #plt.show()	        
    plt.savefig(outfile)
    cp._info('save to %s' % outfile)


if __name__ == '__main__':
        cp.dispatch(__name__)
