# -*- coding: utf-8 -*-
from copy import deepcopy
from itertools import groupby
#from protein import protein

__all__=['msa']

class msa(object):
    def __init__(self, msa_file, protein):
        
        self.name = msa_file
        self.protein = protein
        self.pfam = msa_file[5:12]
        self.msaArray=[]
        self.fa=self.fasta_iter(msa_file)
        for s in self.fa:
            self.msaArray.append(s)
        self.idxMap = self.getIdxMap()
        self.weight = 1.0/len(self.msaArray)
        

            

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
    
    def getWeight(self):
        print len(self.msaArray)
        print self.weight

    def dumpMSA(self):
        print self.name
        for seq in self.msaArray:
            print seq[0]
            print seq[1]+'\n'
            

    # return sequence        
    def getPDBSeq(self):
        s=self.msaArray[0]
        return s[1]
            
            
    def getIdxMap(self):
        # get pdb sequence 
        pdbseq=self.protein.seq
        pdbheader=self.protein.pdb
        
        # get aligned pdb sequence
        mss=self.msaArray[0]        
        msaseq=mss[1]
        msaheader=mss[0]
#        print ph+'\n'+ps
#        print mh+'\n'+ms
        if pdbheader not in msaheader:
            print "getIdxMap()::error:header not match"
            print pdbseq
            print msaseq  
            return 
        
        idxMap={}
        token = -1
        for i in xrange(0,len(pdbseq)):
            for j in xrange(token+1,len(msaseq)):
                if pdbseq[i]==msaseq[j]:
                    idxMap[i]=j
                    token = j
                    break
        if len(idxMap)!=len(pdbseq):
            print "getIdxMap()::error:length not match"
            print len(idxMap)
            print len(pdbseq)
            return
            
        return idxMap            
            
    # return aligned str and aligned indices in string format    
    def getAlignCluster(self, pdbc, msaseq):
        pdbidx = pdbc.getPDBidxArray()
        alignstr = ''
        alignidx = ''
        for i in xrange(0,len(pdbidx)):
            idx = self.idxMap[int(pdbidx[i])]
            alignstr = alignstr + msaseq[idx]
            alignidx = alignidx + ' ' + str(idx)
        return alignstr, alignidx   
        

    def getAlignClusterColumn(self):
        msaclusters=[]
        for i in xrange(0, len(self.protein.clusters)):
#        for i in xrange(0, 2):
            for j in xrange(0, len(self.msaArray)):
                c = deepcopy((self.protein.clusters[i]))
                ret = self.getAlignCluster(c, (self.msaArray[j])[1])
#                print ret
                c.pfam = self.pfam
                c.alignstr = ret[0]
                c.alignidx = ret[1]
                c.seqheader = (self.msaArray[j])[0]
                c.weight = self.weight
#                print c.seqheader
                if j!=0:                
                    c.flag = 1
                msaclusters.append(c)
        return msaclusters

    def writeAlignClusters(self, outfile):
        fp = open(outfile, 'w')
        msaclusters = self.getAlignClusterColumn()
        for c in msaclusters:
#            print c.toString()
            fp.write(c.toString())
        fp.close()

    def printAlignClusters(self):
        msaclusters = self.getAlignClusterColumn()
        print len(msaclusters)
        for c in msaclusters:
            print c.toString()

    def printName(self):
        print '%s: %d %d\n' % (self.name, self.protein.getClusterNum(), len(self.getAlignClusterColumn()))
        
    def writeAlignClustersdb(self, cur):
        msaclusters = self.getAlignClusterColumn()
        for c in msaclusters:
            sql = 'INSERT INTO t_clusters VALUES '+c.toDBString()+';'
#            print "executing: ["+sql+"]"
#            break
            cur.execute(sql) 


