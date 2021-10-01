import sys
import commp as cp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clr
from scipy.cluster.hierarchy import dendrogram
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform
from sklearn import metrics

colorscheme1 = ['#a93a28', '#afc8cd', '#266674', '#fb8c32', '#cbc96d', 
'#60e6c1', '#d7295e', '#008ed0', '#747474']
#1c225c

# 1: tsv data file
# 2: index for x,y,color value
# 3: outfile name
def scatter(args):
    assert len(args) == 5, 'Usage: python utils_vis_sm.py scatter datafile idx.x idx.y idx.color outfile'
    infile = args[0]
    xid = int(args[1])
    yid = int(args[2])
    cid = int(args[3])
    outfile = args[4]

    d = np.loadtxt(infile)
    #plt.scatter(d[xid], d[yid], c=d[cid], cmap="RdYlGn", s=500, edgecolors="black")
    #plt.scatter(d[xid], d[yid], c=d[cid], cmap="RdYlGn")
    #cs = ['#a93a28', '#afc8cd', '#266674', '#fb8c32', '#cbc96d', '#60e6c1', '#d7295e', '#008ed0', '#747474']
    cs = ['#266674','#a93a28']
    mycmap = clr.LinearSegmentedColormap.from_list('mybar', cs, N=256)
    #plt.scatter(d[:,xid], d[:,yid], c=d[:,cid],s=5, facecolors='none', edgecolors='k')
    #plt.scatter(d[:,xid], d[:,yid], c=d[:,cid],s=5, alpha=0.5)
    plt.scatter(d[:,xid], d[:,yid], c=d[:,cid],s=10,alpha=0.8)
    plt.tight_layout()
    #fig.savefig(outfile)
    plt.show()



# input format: 'value','tick_name'
def barplot(args):
    assert len(args) == 2, 'Usage: python utils_vis_sm.py barplot data.vec2 outfile'
    infile = args[0]
    outfile = args[1]

    data = []
    ticks = []

    for line in cp.loadlines(infile):
        sarr = line.split(',')
        data.append(float(sarr[0]))
        ticks.append(sarr[1])

    print(data)
    print(ticks)

    fig, ax = plt.subplots(figsize=(8,6))

    n = len(data)
    index = np.arange(n)
    bar_width = 0.5
    #opacity = 0.4
    outbar = plt.bar(index, data, bar_width, color = colorscheme1[0])
    #plot.xlabel('xlabel')
    #plot.ylabel('ylabel')
    plt.xticks(index + bar_width/2, ticks, rotation='vertical')
    plt.tight_layout()
    fig.savefig(outfile)
    plt.show()



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

    #data = np.loadtxt(infile,delimiter=' ')
    data = np.loadtxt(infile)
    n = len(data)
    #fix, ax = plt.subplots()
    fig, ax = plt.subplots(1,1, figsize=(10,4))
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
    #xlabels = cp.aat01
    xlabels = ['%d' % i for i in index]
    plt.xticks(index + bar_width / 2, xlabels, rotation=90)
    #ax.set_xticklabels(xlabels, rotation=-90)

    plt.legend()

    plt.tight_layout()
    plt.show()
    #plt.savefig(outfile)
    #cp._info('save to %s' % outfile)

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
    ddata = dendrogram(linkage_matrix, labels=xt, leaf_rotation=90, color_threshold=1, link_color_func=lambda x: "k")
    #ddata = _augmented_dendrogram(linkage_matrix, labels=xt, leaf_rotation=90, color_threshold=1)
    #plt.show()	        
    plt.savefig(outfile)
    cp._info('save to %s' % outfile)


def dendrogram_fgdist_1(arglist):
    if len(arglist) < 1:
	    cp._err('Usage: python utils_vis_sm.py dendrogram_fgdist_1 cs3.210dist_1.vec')
    sys.setrecursionlimit(10000)

    datafile = arglist[0]
    outfile = '%s.dendrogram.png' % datafile
    x = np.loadtxt(datafile, delimiter=' ')

    linkage_matrix = linkage(x, "single", metric='correlation')
    plt.figure(figsize=(24, 15))
    #ddata = dendrogram(linkage_matrix)
    xt = cp.aas01
    ddata = dendrogram(linkage_matrix, labels=xt, leaf_rotation=90, color_threshold=1)
    #ddata = _augmented_dendrogram(linkage_matrix, labels=xt, leaf_rotation=90, color_threshold=1)
    #plt.show()	        
    plt.savefig(outfile)
    cp._info('save to %s' % outfile)


def dendrogram_emboss_sm(arglist):
    if len(arglist) < 1:
	    cp._err('Usage: python utils_vis_sm.py dendrogram_emboss_sm b62.sm.dat')
    sys.setrecursionlimit(10000)

    datafile = arglist[0]
    outfile = '%s.dendrogram.png' % datafile
    x = np.loadtxt(datafile, delimiter=' ')

    linkage_matrix = linkage(x, "single", metric='correlation')
    plt.figure(figsize=(24, 15))
    xt = cp.aa201 # emboss sm order
    ddata = dendrogram(linkage_matrix, labels=xt, leaf_rotation=90, color_threshold=1, link_color_func=lambda x: "k")
    #ddata = _augmented_dendrogram(linkage_matrix, labels=xt, leaf_rotation=90, color_threshold=1)
    #plt.show()	        
    plt.savefig(outfile)
    cp._info('save to %s' % outfile)



def dendrogram_mat(arglist):
    if len(arglist) < 3:
        cp._err('Usage: python utils_vis_sm.py dendrogram_mat matfile tickfile outfile')
    matfile = arglist[0]
    tickfile = arglist[1]
    outfile = arglist[2]

    mat = np.loadtxt(matfile, delimiter= ' ')
    xt = cp.loadlines(tickfile)

    dists = squareform(mat)
    #linkage_matrix = linkage(dists, "single")
    linkage_matrix = linkage(dists, "complete")
    #plt.figure(figsize=(25, 28))
    plt.figure(figsize=(15, 15))
    ddata = dendrogram(linkage_matrix, labels=xt, leaf_rotation=90, color_threshold=1, link_color_func=lambda x: "k")
    plt.xticks(fontsize=10)
    plt.tight_layout()
    plt.savefig(outfile)
    #plt.show()

    '''
    dists = squareform(mat)
    linkage_matrix = linkage(dists, "single")
    dendrogram(linkage_matrix, labels=["0", "1", "2"])
    plt.title("test")
    plt.show()
    '''

def dendrogram_linkage(args):
    assert len(args) == 3, 'Usage: python utils_vis_sm.py dendrogram_linkage Z_mat (return from scipy linkage) tick outfile'
    matfile = args[0]
    tickfile = args[1]
    outfile = args[2]
    mat  = np.loadtxt(matfile)
    xt = cp.loadlines(tickfile)
    plt.figure(figsize=(8,8))
    dd = dendrogram(mat, labels=xt, color_threshold=1, link_color_func=lambda x: 'k')  # leaf_rotation=90
    plt.xticks(fontsize=8)
    plt.tight_layout()
    plt.savefig(outfile)
    plt.show()


# plot side-by-side barplot from the data in a two-column file
# total_num: total row number
# group_num: how many barplots to split
def bar_sbs(arglist):
    if len(arglist) < 4:
        cp._err('Usage: python utils_vis_sm.py bar_sbs data.2.vec total_num group_num xtickfile')
    datafile = arglist[0]
    totalnum = int(arglist[1])
    spnum = int(arglist[2])
    n_groups=itv= (totalnum / spnum) +1

    print n_groups

    # x-tick
    xt = cp.loadlines(arglist[3])
    #xt = ['%d' % i for i in range(1,92)]

    data = np.loadtxt(datafile, delimiter=' ')
    d1 = data[:,0]
    d2 = data[:,1]

    fig, ax = plt.subplots(spnum, figsize=(12,8), sharey=True)
    bar_width = 0.25
    opacity = 1.0
    for i in range(0, spnum):
            index = np.arange(n_groups if (i+1)*itv <totalnum else (totalnum - i*itv))
            sxt = xt[i*itv:((i+1)*itv if (i+1)*itv <totalnum else totalnum)]

            ax[i].bar(index, d1[i*itv:((i+1)*itv if (i+1)*itv < totalnum else totalnum)], bar_width,
                            alpha=opacity,
                            color='#75a8b9',
                            label='Ordered region')

            ax[i].bar(index + bar_width, d2[i*itv:((i+1)*itv if (i+1)*itv < totalnum else totalnum)], bar_width,
                            alpha=opacity,
                            color='#8d5543',
                            label='Disordered region')
            '''
            ax[i].bar(index + 2*bar_width, scsc23a[i*itv:((i+1)*itv if (i+1)*itv <210 else 210)], bar_width,
                            alpha=opacity,
                            color='#f3db81',
                            label='SCSC2')
            '''
            plt.sca(ax[i])
            plt.xticks(index + bar_width / 2, sxt, rotation='vertical')

            # plt.xlabel('xxx')
            plt.ylabel('Normalized correlation level')
            plt.legend()

    plt.tight_layout()
    plt.show()

   # fig.savefig(outname)


# generate histogram from data.list
# input: one column (float number) data
def hist(arglist):
    if len(arglist) <2:
        cp._err('Usage: python utils_vis_sm.py hist data.list bin')
    datafile = arglist[0]
    n_bins = int(arglist[1])
    x = np.loadtxt(datafile)
    n = x.shape[0]
    print n
    fig, axs = plt.subplots(tight_layout=True)
    axs.hist(x, bins=n_bins)
    plt.show()

# generate side-by-side histogram from data2.list
# input: two one-column (float number) data files
def hist_sbs(arglist):
    if len(arglist) <2:
        cp._err('Usage: python utils_vis_sm.py hist_sbs data1.list data2.list')
    datafile1 = arglist[0]
    datafile2 = arglist[1]
    data1 = np.loadtxt(datafile1)
    data2 = np.loadtxt(datafile2)
    print data1.shape
    print data2.shape
    x = data1
    y = data2
    n_bins = 50
    fig, axs = plt.subplots(1, 2, sharey=True, tight_layout=True)
    axs[0].hist(x, bins=n_bins)
    axs[1].hist(y, bins=n_bins)
    plt.show()
    #fig.savefig('out.png')
    #cp._info('save to out.png')


# generate simple signal plot (using plot) from two-column (space separted columns) file
# column 1: x-tick
# column 2: values
# originally for visualizing sliding window results of IDP project
def signalplot(arglist):
    if len(arglist) < 2:
        cp._err('Usage: python utils_vis_sm.py singalplot infile')

    infile = arglist[0] 
    outfile = arglist[1]

    valuelist= []
    ticklist = []
    for line in cp.loadlines(infile):
        sarr = line.split(' ')
        ticklist.append(int(sarr[0]))
        valuelist.append(float(sarr[1]))

    print '%d records loaded.' % (len(valuelist))

    fig, ax = plt.subplots(figsize=(14,3),tight_layout=True)
    ax.plot(ticklist, valuelist, label="score")
    #ax.plot(valuelist, label="score")
    ax.grid(True)
    ax.legend()
    plt.axvline(31, linestyle='--', c='r')
    plt.xticks(rotation='vertical', fontname = "monospace")

    #plt.show()
    plt.savefig(outfile)


# batch generate plots for disordered protein
def signalplotws(arglist):
    if len(arglist) < 2:
        cp._err('Usage: python utils_vis_sm.py singalplot_ws infile outfile')

    infile = arglist[0] 
    sarr = infile.split('-')
    dstart = int(sarr[4])
    dend = int(sarr[5])
    outfile = arglist[1]

    valuelist= []
    ticklist = []
    for line in cp.loadlines(infile):
        sarr = line.split(' ')
        ticklist.append(int(sarr[0]))
        valuelist.append(float(sarr[1]))

    #print '%d records loaded.' % (len(valuelist))

    fig, ax = plt.subplots(figsize=(14,3),tight_layout=True)
    ax.plot(ticklist, valuelist, label="score")
    ax.grid(True)
    ax.legend()
    print outfile, dstart, dend
    plt.axvline(dstart, linestyle='--', c='r')
    plt.axvline(dend, linestyle='--', c='r')

    plt.xticks(rotation='vertical', fontname = "monospace")

    plt.show()
    #plt.savefig(outfile)

# plot roc curves 
# infile: true/false score1 score2 ...
def roc(args):
    assert len(args) == 3, 'Usage: python utils_vis_sm.py roc roc_score_file legend,legend2,.. outfile'
    infile = args[0]
    legends = [s for s in args[1].split(',')]
    outfile = args[2]

    data = np.loadtxt(infile)
    print(data.shape)
    p = data.shape[1] # get how many scores

    y = data[:,0]
    curves = [metrics.roc_curve(y, data[:,i], pos_label=1) for i in xrange(1,p)]

    plt.figure(1)
    plt.plot([0, 1], [0, 1], 'k--')
    for k in xrange(0, len(curves)):
        c = curves[k]
        plt.plot(c[0], c[1], color=colorscheme1[k], label=legends[k])

    plt.xlabel('False positive rate')
    plt.ylabel('True positive rate')
    plt.title('ROC curve')
    plt.legend(loc='best')
    plt.savefig(outfile)
    print('save to %s' % outfile)
    plt.show()

def _gridheatmap(data, row_labels, col_labels, ax=None,
            cbar_kw={}, cbarlabel="", **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Arguments:
        data       : A 2D numpy array of shape (N,M)
        row_labels : A list or array of length N with the labels
                     for the rows
        col_labels : A list or array of length M with the labels
                     for the columns
    Optional arguments:
        ax         : A matplotlib.axes.Axes instance to which the heatmap
                     is plotted. If not provided, use current axes or
                     create a new one.
        cbar_kw    : A dictionary with arguments to
                     :meth:`matplotlib.Figure.colorbar`.
        cbarlabel  : The label for the colorbar
    All other arguments are directly passed on to the imshow call.
    """

    if not ax:
        ax = plt.gca()

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    # We want to show all ticks...
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_yticks(np.arange(data.shape[0]))
    # ... and label them with the respective list entries.
    ax.set_xticklabels(col_labels)
    ax.set_yticklabels(row_labels)


    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=True, bottom=False,
                   labeltop=True, labelbottom=False)

    # Rotate the tick labels and set their alignment.
    #plt.setp(ax.get_xticklabels(), rotation=-30, ha="right", rotation_mode="anchor")

    # Turn spines off and create white grid.
    for edge, spine in ax.spines.items():
        spine.set_visible(False)

    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar

def _grid_annotate_heatmap(im, data=None, valfmt="{x:.2f}",
                     textcolors=["black", "white"],
                     threshold=None, **textkw):
    """
    A function to annotate a heatmap.

    Arguments:
        im         : The AxesImage to be labeled.
    Optional arguments:
        data       : Data used to annotate. If None, the image's data is used.
        valfmt     : The format of the annotations inside the heatmap.
                     This should either use the string format method, e.g.
                     "$ {x:.2f}", or be a :class:`matplotlib.ticker.Formatter`.
        textcolors : A list or array of two color specifications. The first is
                     used for values below a threshold, the second for those
                     above.
        threshold  : Value in data units according to which the colors from
                     textcolors are applied. If None (the default) uses the
                     middle of the colormap as separation.

    Further arguments are passed on to the created text labels.
    """

    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = im.norm(data.max())/2.

    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    kw = dict(horizontalalignment="center",
              verticalalignment="center")
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            kw.update(color=textcolors[im.norm(data[i, j]) > threshold])
            text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
            texts.append(text)

    return texts

# heatmap procedure
# colorscheme1 = ['#a93a28', '#afc8cd', '#266674', '#fb8c32', '#cbc96d', '#60e6c1', '#d7295e', '#008ed0', '#747474']
def heatmap(args):
    if len(args) < 3:
        cp._err('Usage: python utils_vis_sm.py heatmap data_npmat.txt xtick.vec ytick.vec range={-1,1}')

    data = np.loadtxt(args[0])
    assert len(data.shape) != 1, 'Error, data has only one row/column'
        
    xticktext = cp.loadlines(args[1]) if args[1]!='na' else np.arange(data.shape[1])
    yticktext = cp.loadlines(args[2]) if args[2]!='na' else np.arange(data.shape[0])

    #minmax = [-0.40, 0.40]
    minmax = [-1.00, 1.00]


    fig, ax = plt.subplots(1,1, figsize=(10,10))

    # customized continuous color bar
    #mycmap = clr.LinearSegmentedColormap.from_list('mybar', ['#266674','#ffffff','#a93a28'], N=256)
    #mycmap = clr.LinearSegmentedColormap.from_list('mybar', ['#ffffff','#a93a28'], N=256)
    mycmap = clr.LinearSegmentedColormap.from_list('mybar', ['#a93a28','#ffffff'], N=256)
    if len(args) == 4:
        minmax = [float(v) for v in args[3].split(',')]
        im = ax.imshow(data, cmap=mycmap, vmin=minmax[0], vmax=minmax[1])
    else:
        im = ax.imshow(data, cmap=mycmap)

    #im = ax.imshow(data, cmap=mycmap, vmin=-1.0, vmax=1.0)
    #im = ax.pcolormesh(data, cmap=mycmap, vmin=-1.0, vmax=1.0)

    # inversed rdbu
    #im = ax.imshow(data, cmap='RdBu_r', vmin=-1.0, vmax=1.0)

    # dirty way but simple way to align the colobar with the plot
    # You can correct for the case where image is too wide using this trick: im_ratio = data.shape[0]/data.shape[1] \
    # plt.colorbar(im,fraction=0.046*im_ratio, pad=0.04) where data is your image
    #cbar = ax.figure.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    # another way to align colorbar with plot, but the plot will no longer be square
    #ax.set_aspect('auto')
    #cbar = ax.figure.colorbar(im, ax=ax)

    # most accurate way to align the colorbar with plot
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = ax.figure.colorbar(im, cax=cax)

    cbar.ax.set_ylabel('colorbar', rotation=-90, va="bottom")

    # show every 10 count of ticks
    ax.set_xticks(np.arange(1, data.shape[1]+1,10))
    ax.set_yticks(np.arange(1, data.shape[0]+1,10)) 
    #ax.set_xticks(np.arange(1, data.shape[1]+1))
    #ax.set_yticks(np.arange(1, data.shape[0]+1)) 
    # ... and label them with the respective list entries.
    ax.set_xticklabels([xticktext[i-1] for i in np.arange(1, data.shape[1]+1,10)], rotation=90)
    ax.set_yticklabels([yticktext[i-1] for i in np.arange(1, data.shape[0]+1,10)])
    #ax.set_xticklabels([xticktext[i-1] for i in np.arange(1, data.shape[1])])
    #ax.set_yticklabels([yticktext[i-1] for i in np.arange(1, data.shape[0])])
    #ax.set_yticklabels(yticktext) 

    # Let the horizontal axes labeling appear on top.

    #plt.xticks(rotation=90)
    ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)    

    fig.tight_layout()
    outfile = args[0]+'.png'
    fig.savefig(outfile)
    cp._info('save figure to %s' % outfile)
    #plt.show()                        


def heatmap_colorgrid(args):
    if len(args) < 3:
        cp._err('Usage: python utils_vis_sm.py heatmap data_npmat.txt xtick.vec ytick.vec range={-1,1}')

    data = np.loadtxt(args[0])
    assert len(data.shape) != 1, 'Error, data has only one row/column'
        
    xticktext = cp.loadlines(args[1]) if args[1]!='na' else np.arange(data.shape[1])
    yticktext = cp.loadlines(args[2]) if args[2]!='na' else np.arange(data.shape[0])

    #minmax = [-1.0, 1.0]

    fig, ax = plt.subplots(1,1, figsize=(10,10))
    mycmap = clr.LinearSegmentedColormap.from_list('mybar', ['#266674','#ffffff','#a93a28'], N=256)
    # two overlap
    mycmap = clr.ListedColormap(['white', '#266674', '#afc8cd', '#a93a28'])
    bounds=[0,0.5,1.0,1.5,2.0]
    norm = clr.BoundaryNorm(bounds, mycmap.N)

    if len(args) == 4:
        minmax = [float(v) for v in args[3].split(',')]
        im = ax.imshow(data, interpolation='nearest', cmap=mycmap, norm=norm, vmin=minmax[0], vmax=minmax[1])
    else:
        im = ax.imshow(data, interpolation='nearest', cmap=mycmap, norm=norm)



    # customized continuous color bar
    #mycmap = clr.LinearSegmentedColormap.from_list('mybar', ['#266674','#ffffff','#a93a28'], N=256)
    #im = ax.imshow(data, cmap=mycmap, vmin=minmax[0], vmax=minmax[1])


    # three overlap
    '''
    mycmap = clr.ListedColormap(['white', '#266674', '#afc8cd', 'red', '#cbc96d', '#008ed0', '#fc9d54', 'cyan'])
    bounds=[0,0.5,1.0,1.5,2.0,2.5,3.0,3.5]
    norm = clr.BoundaryNorm(bounds, mycmap.N)
    '''
    

    #im = ax.imshow(data, interpolation='nearest', cmap=mycmap, norm=norm, vmin=minmax[0], vmax=minmax[1])

    # most accurate way to align the colorbar with plot
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = ax.figure.colorbar(im, cax=cax)

    cbar.ax.set_ylabel('colorbar', rotation=-90, va="bottom")

    # show every 10 count of ticks
    ax.set_xticks(np.arange(1, data.shape[1]+1,10))
    ax.set_yticks(np.arange(1, data.shape[0]+1,10)) 
    #ax.set_xticks(np.arange(1, data.shape[1]+1))
    #ax.set_yticks(np.arange(1, data.shape[0]+1)) 
    # ... and label them with the respective list entries.
    ax.set_xticklabels([xticktext[i-1] for i in np.arange(1, data.shape[1]+1,10)], rotation=90)
    ax.set_yticklabels([yticktext[i-1] for i in np.arange(1, data.shape[0]+1,10)])
    #ax.set_xticklabels([xticktext[i-1] for i in np.arange(1, data.shape[1])])
    #ax.set_yticklabels([yticktext[i-1] for i in np.arange(1, data.shape[0])])
    #ax.set_yticklabels(yticktext) 

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)    

    fig.tight_layout()
    plt.show()                        


def colorpalette(args):
    #      red        lightblue  darkgreen  orange    dullyellow  lightgreen magenta   lightblue  gray
    cs = ['#a93a28', '#afc8cd', '#266674', '#fb8c32', '#cbc96d', '#60e6c1', '#d7295e', '#008ed0', '#747474']
    # start showing colorbar
    #cs = ['#266674','#ffffff','#a93a28']

    # plot psudodata
    csdata = np.array([range(len(cs))])
    fig, ax = plt.subplots(figsize=(8,5))
    #plt.gca().set_visible(False) # hide it
    plt.gca().get_yaxis().set_visible(False)

    ax.set_xticks(np.arange(0, csdata.shape[1]+1))
    ax.set_yticks(np.arange(0, csdata.shape[0]+1)) 
    # ... and label them with the respective list entries.
    ax.set_xticklabels(cs)
    #ax.set_yticklabels([yticktext[i-1] for i in np.arange(1, data.shape[0]+1,10)])

    mycmap = clr.LinearSegmentedColormap.from_list('mybar', cs, N=256)

    im = ax.imshow(psudodata, cmap=mycmap)
    cbar = ax.figure.colorbar(im, ticks=[range(len(cs))], orientation="horizontal")
    cbar.ax.set_xticklabels(cs, rotation=-90)

    plt.show()
    #plt.savefig("colorbar.pdf")    

class seq(object):
    def __init__(self, props):
        self.resn = props['resn']
        self.resi = props['resi']
        self.ss = props['ss']
        self.highlight = props['highlight']
        self.consv = props['consv']

    def __repr__(self):
        return ('%s %s %s %s') % \
            (
                self.resn,
                self.resi,
                self.highlight,
                self.consv
            )

def t_seq(args):
    res = {'resi': '30', 'resn': 'LYS' }
    print res['resi']
    sq = seq(res)
    print repr(sq)

def seqdraw(args):
    # import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6,50))

    def __parseseq(line):
        sarr = line.split(' ')
        return {
            'resn': sarr[0], 'resi': sarr[1], 'ss': sarr[2], 'highlight': int(sarr[3]), 'consv': sarr[4]
            }

    #seqs = [seq(__parseseq(line)) for line in cp.loadlines('t_seqdraw.data')]
    seqs = [seq(__parseseq(line)) for line in cp.loadlines('3k54_seqdraw.data')]
    print len(seqs)

    sscolor = {'s':'#ff39ff', 'h':'#39ffff', 'l': '#ff9999'}

    fs = 8 

    # draw seq
    for i in range(len(seqs)):
        x = 0 
        y = float(len(seqs)-i-1)/len(seqs)
        sq = seqs[i]
        # draw sequence
        ax.text( 
            x,y, 
            '%s %s' % (sq.resn, sq.resi), 
            ha='center', va='center',
            fontsize=fs, family= 'monospace', color='black', 
            bbox=dict(pad=0.3, facecolor=sscolor[sq.ss], edgecolor=sscolor[sq.ss], alpha=0.3, boxstyle='square'))
        # highlight
        x+=0.07
        fc = 'gray' if sq.highlight == 1 else 'white'
        ax.text( 
            x,y, 
            ' ',
            ha='center', va='center',
            fontsize=fs-2, family= 'monospace', color='black', 
            bbox=dict(pad=0.3, facecolor=fc, edgecolor='white', boxstyle='square'))
        # conservation
        x+=0.05
        ax.text( 
            x,y, 
            '%s' % (sq.consv) if sq.consv != 'na' else ' ',
            ha='center', va='center',
            fontsize=fs-1, family= 'monospace', color='black', 
            bbox=dict(pad=0.3, facecolor='white', edgecolor='gray', alpha=0.3, boxstyle='square'))
    plt.axis('off')
    #plt.show()
    fig.tight_layout()
    plt.savefig("3k54.png")

    '''
    from matplotlib.offsetbox import (TextArea, DrawingArea, OffsetImage, AnnotationBbox)
    #xy = (0.5, 0.7)
    xy = (0, 0)
    ab = AnnotationBbox(TextArea("textarea"), xy,
                    xybox=([0, 0]), # text location
                    xycoords='data',
                    boxcoords=("axes fraction", "data"),
                    box_alignment=(0., 0.5),
                    )
    ax.add_artist(ab)
    #xy = (0, 0.5)
    ab = AnnotationBbox(TextArea("textarea1"), xy,
                    xybox=([0, 0.5]),
                    xycoords='data',
                    boxcoords=("axes fraction", "data"),
                    box_alignment=(0., 0.5),
                    )
    ax.add_artist(ab)
    '''
    '''
    ax.text(
        0,1, # text location, 0:bottom, 1:top
        'foobar', # text
        fontsize=12, # fontsize # equal to "size=12"
        color='black', # fontcolor
        bbox=dict(facecolor='red', edgecolor='red', boxstyle='round') # box attributes
        )
    t = ax.text(
        0, 0, "Direction", ha="center", va="center", rotation=45, size=15,
        bbox=dict(boxstyle="rarrow,pad=0.3", fc="cyan", ec="b", lw=2))
    '''

if __name__ == '__main__':
        cp.dispatch(__name__)
