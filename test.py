#!/usr/bin/python
from protein import protein
from cluster import cluster
from msa import msa
from AAShingle import AAShingle
from sklearn.cluster import spectral_clustering
import sys

def main():
	graph=''
	labels = spectral_clustering(graph, n_clusters=4)
####################################### test Tips
#    p=protein('1x9d.tip', 'alpha',center='TIP')
    #print len(p.atoms)
#    p.filterClusters()
    #print len(p.clusters)
   # print p.center
#    p.printPDB()
   # p.writeChainATips('AAtips.txt', '1qgr.tip')








    
####################################### test AAShingle
#    shingle = AAShingle()
#    print shingle.getValue('A')
#    print shingle.shingle2Index('RAI')
#    print shingle.index2Shingle(5607)
#    print shingle.shingle2Index('RAV')
#    print shingle.index2Shingle(5617)    
#    
#    
#    s='IKYILGKISAMRKEMCEKYEKCENSKEVLAENNLNLPKMAEKDGCFQSGFNQETCLMRITTGLVEFQIYLDYLQKEYESNKGNVEAVQISTKALIQTLRQKGKNPDKATTPNPTTNAGLLDKLQSQNEWMKNTKIILILRSLEDFLQFSLRAIR'
##    s='IKYILGKI'
#    print len(s)
#    print s
#    print '['+shingle.seq2Shingle(s,3)+']'
#    
#    shingle.shingleStr2Fasta('PF00489_seed.txt','PF00489_seed.shingle',3)
#    shingle.shingleIndex2Fasta('PF00489_seed.txt','PF00489_seed.index',3)
####################################### test msa    
    
####################################### test msa
#    msafile = sys.argv[1]
#    pdbfile = (sys.argv[1])[0:4]+'.pdb'
#    topology = sys.argv[2]
#    
#    print msafile
#    print pdbfile
#    print topology
##    p=protein('1alu.pdb', 'alpha')
#    p=protein(pdbfile, topology)
#    p.filterByCA()
##    print p.getSeq()
#
##    print p.clusters[0].dump()
#    
##    m=msa('1alu_PF00489.fa', p)
#    m=msa(msafile, p)
#    m.printAlignClusters()
#    m.printName()
#########################################    
#    m.writeAlignClusters('1alu_PF00489.clusters')    
#    m.getWeight()
#    print m.pfam
#    print (m.msaArray[0])[1]
#    print '\n'
    #test msa index mapping

    
    

#    p.dumpClusters()
#    msac = m.getAlignClusterColumn()
#    print len(p.clusters)
#    print len(msac)

#    print p.clusters[0].getPDBidxArray()
#    print p.clusters[0].str
#    print (p.clusters[0]).str, (p.clusters[0]).pdbidx
#    print m.getAlignCluster(p.clusters[0], (m.msaArray[0])[1]) 
    
#    msac[0].dump()
#    msac[1].dump()
#    msac[2].dump()
#    msac[3].dump()
#    for c in msac:
#        print c.toString()

#    idxmap = m.getIdxMap()
#    print idxmap
#    
#    print p.seq[153]
#    print ((m.msaArray[0])[1])[249]

#    m.writeAlignClusters('1alu_PF00489.clusters')
#    m.writeAlignClustersdb('1alu_PF00489.clusters')    
#
#    s1=''
#    s2=''
#    
#    pdbseq=p.getSeq()
#    msaseq=m.getPDBSeq()
#
#    for key in idxmap.keys():
#        s1=s1+msaseq[idxmap[key]]
#        s2=s2+pdbseq[key]
#        
#    print s1
#    print s2
#    print pdbseq
#    print msaseq    
#    print idxmap
#    
#    print pdbseq[0]+"-"+msaseq[30]
##    m.dumpMSA()
##
#    s1=m.getPDBSeq()
#    print '>%s\n%s\n' % (s1[0], s1[1])
    # test protein
#    pdb_filename='1l2y.pdb'
#    # pdbname, top='', pfam='', center='CA', cutoff=7, scutoff=2, flag=0, desc='', nbcutoff=5
#    p=protein(pdb_filename, 'cluster')
#    p.getPairwiseOf(19)
##    p.dump()
#    p.filterByCA()
##    p.dumpClusters()
#    p.clusters[0].dump()
#    print p.clusters[0].getPDBidxArray()

#    pdb_filename='1l2y.pdb'
#    p=protein(pdb_filename,'test')
##    p.writeSeq('1l2y.fa')
##    p.writeChainACA('1alu_ca.pdb')
##    
#
#    p.pairwise()
#    pd=p.pairwiseDict
#    for k in pd.keys():
#        print "%s: %f" % (k,pd[k])
#    print pd['0-0']       
#    
#    print "%s: %f" % ('0-1',pd['0-1'])
#    print "%s: %f" % ('1-0',pd['1-0'])
#    print len(pd)

    # test cluster
#    c=cluster('1l2y','d1l2ya_', 'pfam_1l2y', 'AABBCC', '1 2 3 4 5 6', 'header_1l2y', 'aabbcc', '11 12 13 14 15 16', 'CA', 9, 2, 0, 'Trp-Cage Miniprotein Construct TC5b')
#    c.dump()
#    s=c.getString()
#    print s

	pass
if __name__=="__main__":
	main()
