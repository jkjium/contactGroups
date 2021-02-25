import sys
import numpy as np

# parse sm (substituion matrix) file
class smatrix(object):
    # read emboss format matrix or 20x20 standard matrix from a file
    def __init__(self, smfile):
        self.name = smfile
        self.aa = []
        scorelist = []
        with open(smfile) as fp:
            for line in fp:
                line = line.strip()
                if len(line)== 0 or line[0] == '#':
                    continue
                if any(c.isdigit() for c in line): # score line
                    smline = line.split()
                    scorelist.append([int(i) for i in smline[1:]])
                else: # alphabet list
                    self.aa = line.split()

        npscore = np.array(scorelist)
        self.core = npscore[:20,:20] # standard 20 x 20 elements
        self.edge = np.copy(npscore) # B Z X *
        self.score =dict(('%s%s' % (self.aa[i],self.aa[j]), self.core[i][j]) for i in range(20) for j in range(20))
    
def sm2mafft():
    alphabets = {
        'A': '0x41', 'C': '0x43', 'D': '0x44', 'E': '0x45', 'F': '0x46', 
        'G': '0x47', 'H': '0x48', 'I': '0x49', 'K': '0x4b', 'L': '0x4c',
        'M': '0x4d', 'N': '0x4e', 'P': '0x50', 'Q': '0x51', 'R': '0x52',
        'S': '0x53', 'T': '0x54', 'V': '0x56', 'W': '0x57', 'Y': '0x59'        
    }
    smfile = sys.argv[1]
    sm =smatrix(smfile)
    for i in range(20):
        for j in range(i,20):
            print("%s %s %d" % (alphabets[sm.aa[i]], alphabets[sm.aa[j]], sm.core[i][j]))

if __name__=='__main__':
    sm2mafft()