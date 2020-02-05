import commp as cp
from scipy.stats import skew

def foo(arglist):
    print arglist

# input: column infile, column index(0-based), outfile
# output: .vec file with unscrewed values (+zscore)
# use cubic root to unskew left screwd distribution
def zscoredist(arglist):
    if len(arglist) < 4:
        cp._err('Usage: python utils_flaten.py zscoredist PF00240.rdca 8 PF00240.trans.dca.vec 1')
    infile = arglist[0]
    idx = int(arglist[1]) # column index
    outfile = arglist[2]
    transflag = int(arglist[3])

    rawlist = []
    for line in cp.loadlines(infile):
        sarr = line.split(' ')
        rawlist.append(float(sarr[idx]))
    oldskew = skew(rawlist)
    translist = [x ** ( 1. / 3 ) for x in dcalist]
    newskew = skew(translist)

if __name__=='__main__':
    cp.dispatch(__name__)