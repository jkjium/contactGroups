# -*- coding: utf-8 -*-
from itertools import groupby

__all__=['AAShingle']

class AAShingle(object):# -*- coding: utf-8 -*-
#sort aa.list | awk '{print "self.baseAlphabet['\''"$1"'\'']="NR}'

    def __init__(self):
        self.baseAlphabet={}
        self.baseAlphabet['A']=0
        self.baseAlphabet['C']=1
        self.baseAlphabet['D']=2
        self.baseAlphabet['E']=3
        self.baseAlphabet['F']=4
        self.baseAlphabet['G']=5
        self.baseAlphabet['H']=6
        self.baseAlphabet['I']=7
        self.baseAlphabet['K']=8
        self.baseAlphabet['L']=9
        self.baseAlphabet['M']=10
        self.baseAlphabet['N']=11
        self.baseAlphabet['P']=12
        self.baseAlphabet['Q']=13
        self.baseAlphabet['R']=14
        self.baseAlphabet['S']=15
        self.baseAlphabet['T']=16
        self.baseAlphabet['V']=17
        self.baseAlphabet['W']=18
        self.baseAlphabet['Y']=19

    def getValue(self, A):
        return self.baseAlphabet[A]
    
    
    def shingle2Index(self, s, base=20):
        value = 0
        n = len(s) - 1
        for i in s:         
               value = value  + self.baseAlphabet[i]*pow(base,n)
               n = n -1 
        return value
 
 
    def index2Shingle(self, index, base=20):
        return ((index == 0) and  "0" ) or ( self.index2Shingle(index // base, base).lstrip("0") + "ACDEFGHIKLMNPQRSTVWY"[index % base])
#       shingle = ""
#       while index != 0:
#           remainder = index % base
#           print remainder
#           if 20 > remainder:
#               remainder_string = self.baseAlphabet[remainder]
#           shingle = remainder_string+shingle
#           index = index / base
#       return shingle    
    def fasta_iter(self, fasta_name):
        """
        given a fasta file. yield tuples of header, sequence
        """
        fh = open(fasta_name)
        # ditch the boolean (x[0]) and just keep the header or sequence since
        # we know they alternate.
        faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
        for header in faiter:

            # drop the ">"
            header = header.next()[1:].strip()
            # join all sequence lines to one.
            seq = "".join(s.strip() for s in faiter.next())
            yield header, seq
            
            
    # s: sequence
    # n: shingle length 
    def seq2Shingle(self, s, n):
        singles=''
        for i in xrange(0, len(s)-n+1):
            singles = singles +' '+ s[i:i+n]
        return singles.lstrip(' ')
            
            
    # infile: fasta file (not alignment)
    # outfile: output file with fasta foramt
    # n: shingle length
    def shingleStr2Fasta(self, infile, outfile, n):
        fa=self.fasta_iter(infile)
        fh=open(outfile, 'w')
        for s in fa:
            fh.write('>'+s[0]+'\n')
            fh.write(self.seq2Shingle(s[1],n)+'\n')
        fh.close()

        
    # infile: fasta file (not alignment)
    # outfile: output file with fasta foramt
    # n: shingle length
    def shingleIndex2Fasta(self, infile, outfile, n):
        fa=self.fasta_iter(infile)
        fh=open(outfile, 'w')
        for s in fa:
            fh.write('>'+s[0]+'\n')
            vstr = ''
            for sh in self.seq2Shingle(s[1],n).split(' '):
                vstr=vstr+' '+str(self.shingle2Index(sh))
            fh.write(vstr.lstrip(' ')+'\n')
        fh.close()        
        
        
        