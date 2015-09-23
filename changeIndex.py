import sys
def main():
    if len(sys.argv)<3:
        print "Usage changeIndex.py db_export pdbfile"
        return
    dbout = sys.argv[1]
    pdbfile = sys.argv[2]
    #fin = open('1x9d_short.csv', 'r')
    fin = open(dbout, 'r')
    lines = fin.readlines()
    fin.close()

    #fpdb = open('1x9d.pdb' ,'r')
    fpdb = open(pdbfile ,'r')
    rl = fpdb.readlines()	
    fpdb.close()
 
    for i in xrange(0,len(lines)):  
        line = lines[i].strip()
	strArr = line.split(',')
	hamming= sum([ch1!=ch2 for ch1, ch2 in zip(strArr[0],strArr[1])])	
	p=strArr[2].replace('"','').strip().split(' ')
	
	resSeqStr=''
	for j in xrange(0,len(p)):
		atom = rl[int(p[j])]
		resSeq = atom[22:26].strip()
		resSeqStr=resSeqStr+' '+resSeq
   	#print resSeqStr
	print '%d,%s,%s' % (hamming, line, resSeqStr.strip())

    pass
if __name__=="__main__":
	main()# -*- coding: utf-8 -*-


