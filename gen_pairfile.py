import sys
import random

def readseq(seqfile):
	with open(seqfile) as f:
		line = f.readline()
		return line.strip()

# iterate and write all the combination files
def iteratepair():
        if len(sys.argv) < 3:
                print 'Usage: python gen_pairfile.py seqfilelist'
                return

        fin = open(sys.argv[2])
        lines = fin.readlines()
        fin.close()

        count=0
        for i in xrange(0, len(lines)):
                for j in xrange(i+1, len(lines)):
                        seqf1 = lines[i].strip()
                        seqf2 = lines[j].strip()
                        seq1 = readseq(seqf1)
                        seq2 = readseq(seqf2)
                        if len(seq1)==0 or len(seq2)==0:
                                print 'error seq: \n%s\n%s' % (seq1, seq2)
                                exit(1)
                        outfile = 'p.%s-%s.seq' % (seqf1[0:4], seqf2[0:4])
                        fout=open(outfile, 'w')
                        fout.write('%s %s' % (seq1, seq2))
                        fout.close()
                        count+=1
        print '%d files generated.' % count


# output a list of pair files
# without inlcuding any pair from exclude list
# pair name in 'p.%s-%s.seq' format
# exclude line in '%s %s' format
def randompair():
        if len(sys.argv) < 6:
                print 'Usage: python gen_pairfile.py randompair pairfile excludelist number outfile'
                return

	pairfile = sys.argv[2]
        excludefile = sys.argv[3]
	num = int(sys.argv[4])
	outfile = sys.argv[5]

	pairlist = []
	with open(pairfile) as f:
		for line in f:
                        if len(line)<2:
                                continue
			pairlist.append(line.strip())
	print '%d pair names loaded. ' % len(pairlist)

        excludelist = []
        with open(excludefile) as f:
                for line in f:
                        if len(line)<2:
                                continue
                        line = line.strip()
                        nameArray = line.split(' ')
                        name = 'p.%s-%s.seq' % (nameArray[0], nameArray[1])
                        excludelist.append(name)
        print '%d exclude pair names loaded. ' % len(excludelist)

        count=0
	indexlist = random.sample(range(len(pairlist)), num+len(excludelist))
	fout = open(outfile,'w')
	for i in indexlist:
                if pairlist[i] in excludelist:
                        print 'skip pair: %s' % pairlist[i]
                        continue
		fout.write('%s\n' % pairlist[i])
                count+=1
	fout.close()
	
	print '%d pairs written.' % count
		



def main():
        dispatch = {
                'iteratepair':iteratepair, 'randompair':randompair
        }

        if len(sys.argv)<2:
                for k in dispatch:
                        dispatch[k]()
                return

        cmd = sys.argv[1]

        flag = False
        for key in dispatch:
                if key == cmd:
                        dispatch[key]()
                        flag = True
        if flag == False:
                print 'Wrong cmd string: %s' % sys.argv[1]


if __name__ == '__main__':
	main()
