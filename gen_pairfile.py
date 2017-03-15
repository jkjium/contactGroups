import sys
import random

def readseq(seqfile):
	with open(seqfile) as f:
		line = f.readline()
		return line.strip()

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


def randompair():
        if len(sys.argv) < 5:
                print 'Usage: python gen_pairfile.py randompair pairfile number outfile'
                return

	pairfile = sys.argv[2]
	num = int(sys.argv[3])
	outfile = sys.argv[4]

	pairlist = []
	with open(pairfile) as f:
		for line in f:
			pairlist.append(line.strip())
	print '%d pair names loaded. ' % len(pairlist)

	indexlist = random.sample(range(len(pairlist)), num)
	fout = open(outfile,'w')
	for i in indexlist:
		fout.write('%s\n' % pairlist[i])
	fout.close()
	
	print '%d pairs written.' % len(indexlist)
		



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
