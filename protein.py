#!/usr/bin/python
import sys
import math
import copy
import numpy as np

from atom import atom
from AAmap import AAmap
from cluster import cluster
from ncg import ncg
import commp as cp

__all__=['protein']

class protein(object):
    
    
    # read pdb from file
    def __init__(self, pdbname, chain = 'all', top='', pfam='', center='CA', cutoff=5, scutoff=1, flag=0, desc='', nbcutoff=4):
        self.atoms=[]
        #dictionary for pairwise distance
        self.pairwiseDict={}
        self.clusters=[]

        #pdb, top, pfam, str, pdbidx, seqheader, alignstr, alignidx, center, cutoff, scutoff, flag, desc
        #self.pdb = pdbname[len(pdbname)-8:len(pdbname)-4]
        self.pdbfile = pdbname
        self.pdb = pdbname[:-4]
        self.chain = chain
        self.top = top
        self.pfam = pfam
        self.center = center
        self.cutoff = cutoff
        self.scutoff = scutoff
        self.seqheader = self.pdb
        self.flag = flag
        self.desc = desc
        
        self.nbcutoff = nbcutoff
        self.ca = []
        self.reslist = []
        
        fin=open(pdbname, 'r')
        lines=fin.readlines()
        fin.close()
       
        lastname =''
        lastres = ''
        ratoms = [] # store residue atoms
        rlist = []
        aamap = AAmap()
        lastresi = '' # redundant ... should make this better (when have time)
        for i in xrange(0,len(lines)):
            line = lines[i]
            # load only one model
            if 'END' in line[0:6]:
                break
            if line[17:20].strip() not in aamap.AAA2A:
                #cp._info('skipped: non-AA residue name: \n%s' % line)
                continue
            if self.chain != 'all':
                if (self.chain != line[21]):
                    continue
            if line[0:6]=='ATOM  ':
                at = atom(lines[i])
                # at.name = {CA, N, ...}
                if (at.name == lastname) and (at.resSeq == lastres): # discard alternative locations (same residue number and same atom name)
                    #print '[%s]::alter loc:\n%s' % (self.pdbfile, lines[i])
                    #if (line[16]==' ' or line[16]=='A'): # to avoid alternative location
                    continue
                else: # for all the normal ATOM entries
                    self.atoms.append(at)
                    if at.name.strip()=='CA':
                        self.ca.append(at)
                    lastname = at.name
                    lastres = at.resSeq

                # store atoms by residueID 
                if (at.resSeq != lastresi): # reach a new residue
                    if(len(ratoms)!=0):
                        self.reslist.append(ratoms) # put atoms for last residue in reslist
                    ratoms = [] # init a new list to store the current residue atoms
                    lastresi = at.resSeq
                    ratoms.append(at)
                else: # still for the same residue
                    ratoms.append(at)
        # store the last residue
        if(len(ratoms)!=0):
            self.reslist.append(ratoms)

        # map for Chain+Resi : (index in sequence, ResName)
        # 'B529': (132, 'V')
        self.resDict = {} # assigned in self.getSeq() function
        # resAtoms, a list of lists, each (element) list contains atoms of residues
        # resArray, gives a list of keys eg. (A,Q,70), (A,I,71), (A,V,72)
        self.seq, self.resArray, self.resAtoms = self.getSeq()

        # resAtomsDict[key] = list() of atom
        # 1gzh_A.pdb resAtomsDict['A196'] = 
        '''
            ATOM    745  N   ARG A 196      20.463  12.417  35.872  1.00 27.01           N
            ATOM    746  CA  ARG A 196      20.781  11.116  36.427  1.00 29.11           C
            ATOM    747  C   ARG A 196      19.585  10.190  36.233  1.00 31.25           C
            ATOM    748  O   ARG A 196      18.502  10.622  35.803  1.00 31.36           O
            ATOM    749  CB  ARG A 196      21.060  11.231  37.928  1.00 28.76           C
            ATOM    750  CG  ARG A 196      22.275  12.054  38.291  1.00 30.62           C
            ATOM    751  CD  ARG A 196      22.548  11.982  39.782  1.00 29.49           C
            ATOM    752  NE  ARG A 196      23.663  12.845  40.163  1.00 34.64           N
            ATOM    753  CZ  ARG A 196      24.097  12.997  41.412  1.00 36.70           C
            ATOM    754  NH1 ARG A 196      25.120  13.806  41.674  1.00 36.31           N1+
            ATOM    755  NH2 ARG A 196      23.508  12.337  42.403  1.00 38.20           N
        '''
        self.resAtomsDict = {} 
        for rlist in self.resAtoms:
            #key = '%s%d%s' % (rlist[0].chainID, rlist[0].resSeq, rlist[0].iCode)
            # insertion of residues are skipped in previous self.getSeq()

            key = '%s%d' % (rlist[0].chainID, rlist[0].resSeq)
            if key in self.resAtomsDict:
                cp._info('duplicate residue found: %s in %s' % (key, self.pdbfile))
                continue
            self.resAtomsDict[key] = rlist

        # some residue does not have CA!! 1e6i.aln.pdb the last residue
        #aamap = AAmap()
        #self.seq = ''.join([aamap.getAAmap(a.resName) for a in self.ca])

        # map for sequence index: Chain+Resi(ResName)
        # 132 : 'B529(V)'
        self.seqDict = {-1: '.'}
        for r in self.resDict:
            self.seqDict[self.resDict[r][0]] = '%s(%s)' % (r, self.resDict[r][1])
        

    # get atom by atom name, ie. CA, CB, ...
    def atomsbyname(self, aname):
        return [a for a in self.atoms if a.name.strip() == aname]

    # get geometrical center for each residue
    def atomsbygmcenter(self):
        ats = []
        for al in self.resAtoms:
            x=0.0
            y=0.0
            z=0.0

            # save CA as a template
            # get accumulative coordinates
            for a in al:
                x+=a.x
                y+=a.y
                z+=a.z

            # replace geom center to template coordinate 
            # and save for output
            reta = copy.copy(al[0])
            reta.x = x/len(al)
            reta.y = y/len(al)
            reta.z = z/len(al)
            ats.append(reta)

        return ats


    # get geometrical center for each residue side chain
    # except for gly
    def atomsbyscgmcenter(self):
        bb = ['N', 'CA', 'C', 'O']
        ats = []
        for al in self.resAtoms:
            x=0.0
            y=0.0
            z=0.0

            # save CA as a template
            # get accumulative coordinates
            count = 0
            if al[0].resName == 'GLY':
                for a in al:
                    if a.name.strip() == 'CA':
                        count+=1
                        x=a.x
                        y=a.y
                        z=a.z
            else:
                for a in al:
                    #if a.name.strip() in bb and a.resName!='GLY':
                    if a.name.strip() in bb:
                        continue
                    count+=1
                    x+=a.x
                    y+=a.y
                    z+=a.z

            # replace geom center to template coordinate 
            # and save for output
            reta = copy.copy(al[0])
            if count == 0:
                cp._info('err: residue: %s %d %s has no side chain' % (self.pdbfile, reta.resSeq, reta.resName))
                continue
            reta.x = x/count
            reta.y = y/count
            reta.z = z/count
            ats.append(reta)

        return ats        
    
    # return a redundant list of Chain+resi
    # the redundancy is the contact affinity
    def contactbyallatom(self, chain, resi, cutoff, seqcutoff=0.0):
        key = '%s%d' % (chain, resi)
        neighbors = []
        for r in self.resAtomsDict:
            # do not compare with itself
            if key == r:
                continue
            for a in self.resAtomsDict[r]:
                for b in self.resAtomsDict[key]:
                    dist =  np.linalg.norm(np.array((a.x, a.y, a.z))-np.array((b.x, b.y, b.z)))
                    if dist <= cutoff and abs(a.resSeq-b.resSeq) > seqcutoff:
                        neighbors.append('%s%d' % (a.chainID, a.resSeq))
        return neighbors

    # return a redundant list of Chain+resi
    # the redundancy is the contact affinity
    def neighborsbyca(self, chain, resi, cutoff, seqcutoff=0.0):
        key = '%s%d' % (chain, resi)
        neighbors = []
        # locate the target ca atom
        targetid = -1
        for i in range(len(self.ca)):
            a = self.ca[i]
            if a.resSeq == int(resi) and (a.chainID == chain):
                targetid = i
                break
        if targetid == -1:
            cp._err('target resi %d not found in %d' % (int(resi), self.pdbfile))
        b = self.ca[targetid]
        for a in self.ca:
            dist =  np.linalg.norm(np.array((a.x, a.y, a.z))-np.array((b.x, b.y, b.z)))
            if dist <= cutoff and abs(a.resSeq-b.resSeq) > seqcutoff:
                neighbors.append('%s%d' % (a.chainID, a.resSeq))
        return neighbors

    # output all pairwise residue distance
    # distance evalued by any atom
    # the minimum distance is recorded in the result
    # output format:  a list of triplet tuples
    # [(resiA, resiB, minDist), ( .. ) ... ]
    def residistbyallatom(self):
        retlist = []
        # keys in resAtomsDict: 'A100', ...
        rlist = self.resAtomsDict.keys()
        for i in xrange(0, len(rlist)):
            for j in xrange(i+1, len(rlist)):
                # get asymmetrical keys
                r1 = rlist[i]
                r2 = rlist[j]
                # calculate minimum distance among all atoms
                mindist = min([np.linalg.norm(np.array((a.x, a.y, a.z))-np.array((b.x, b.y, b.z))) for a in self.resAtomsDict[r1] for b in self.resAtomsDict[r2]])
                # extract resi number
                resi1 = int(r1[1:])
                resi2 = int(r2[1:])
                # make sure residue with smaller ID is in front in the output
                if resi1 < resi2:
                    retlist.append((r1, r2, mindist))
                else:
                    retlist.append((r2, r1, mindist))
        return retlist



    # output residue contact by dist cutoff and seqcutoff
    # no redundancy 
    def contactbycutoff(self, atomset, cutoff, seqcutoff=0.0):
        cgs = []
        for i in xrange(0, len(atomset)):
            for j in xrange(i+1, len(atomset)):
                a = atomset[i]
                b = atomset[j]
                dist =  np.linalg.norm(np.array((a.x, a.y, a.z))-np.array((b.x, b.y, b.z)))                
                if dist <= cutoff and abs(a.resSeq-b.resSeq) > seqcutoff:
                    cgs.append((a,b))
        return cgs

    # output residue contact by nearest neighbor
    # redundant contact removed
    def contactbynearest(self, atomset, size):
        cgs = []
        ncgArray = []
        for a in atomset:
            c = ncg(a, size)
            ncgArray.append(c)

        # grow by nearest neighbor
        # remove duplicate
        dup = set()
        for c in ncgArray:
            c.grow(atomset)
            key = ' '.join(sorted(['%s%d' % (a.chainID, a.resSeq) for a in c.atoms]))
            #if key in dup:
            #    print 'duplicate key: %s' % key
            if key not in dup:            
                cgs.append(c.atoms)
                dup.add(key)
        return cgs

    
    # print PDB content
    def printPDB(self):
        for i in xrange(0,len(self.atoms)):
            a=self.atoms[i]
            a.dump()
        
    # return a list of 'str' Coordinates
    def printcoor(self):
        return ['%d %.3f %.3f %.3f' % (a.resSeq, a.x, a.y, a.z) for a in self.atoms] 

    # return sequence extracted from pdb file
    # assign values for self.resDict['B641'] = (seqpos, 'R')
    def getSeq(self):
        aamap = AAmap()
        seq=''
        #last_resSeq = -1 # 1a8v the first resi starts from -1 !!!!
        last_resSeq = -9999 # 1a8v the first resi starts from -1 !!!!
        seqPos = 0
        resArray = []

        resAtomsAll = []
        resatoms = []
        insertion_residues={} # record insertion residues
        for i in xrange(0,len(self.atoms)):
            a=self.atoms[i]
            #print '[%s]' % a.iCode
            if last_resSeq != a.resSeq: # this will naturally avoid loading insertion residues, since insertion residues have the same resSeq
                seq=seq+aamap.getAAmap(a.resName)
                last_resSeq = a.resSeq

                key = '%s%d' % (a.chainID, a.resSeq)
                self.resDict[key] = (seqPos, seq[seqPos])
                seqPos+=1

                #resArray.append('%s %s %s' % (a.chainID,aamap.getAAmap(a.resName),str(a.resSeq)))
                resArray.append((a.chainID,aamap.getAAmap(a.resName),a.resSeq))

                if len(resatoms)>0:
                    resAtomsAll.append(resatoms)
                    resatoms=[]

            if a.iCode!=' ':
                cp._info('insertion residue %d%s found in %s' % (a.resSeq, a.iCode, self.pdbfile))
                if a.resSeq not in insertion_residues: 
                    insertion_residues[a.resSeq] = a.iCode
                    resatoms.append(a) 
                else:
                    if insertion_residues[a.resSeq] == a.iCode: # atoms in the same residue
                        resatoms.append(a) 
                    else: # inserted residues of the same resi, discard atom
                        continue
            else:
                resatoms.append(a) # add normal atoms

        # after loop add the last res into resatoms
        # only resSeq change trigger adding above
        if len(resatoms)>0:
            resAtomsAll.append(resatoms)

        return seq, resArray, resAtomsAll       
           
    # print PDB sequence into fasta file
    def writeFASTA(self):
        fafile = self.pdb+'.fa'
        aamap = AAmap()

        seq=''
        count = 0
        last_resSeq = -1
        for i in xrange(0,len(self.atoms)):
            a=self.atoms[i]
            if last_resSeq != a.resSeq:
                seq=seq+aamap.getAAmap(a.resName)
                last_resSeq = a.resSeq
                count+=1
        seq=seq+'\n'
        header = '>%s/1-%d\n' % (self.pdb, count)
        print header+seq

        fp=open(fafile, 'w')
        fp.write(header+seq)
        fp.close()


    # extract all the CA atoms from the pdb file
    def writeCA(self, filename, chain='all'):
        fd=open(filename, 'w')
        count=0
        for i in xrange(0,len(self.atoms)):
            a=self.atoms[i]
            #a.dump()
            if chain == 'all':
                if 'CA' in a.name:
                    fd.write(a.writeAtom())
                    count=count+1
            else:
                if 'CA' in a.name and a.chainID == chain:
                    fd.write(a.writeAtom())
                    count=count+1
        fd.close()
        if count==0:
            print "No atom written in [%s]!" % (filename)      
        
        
    
    # extract all the CA atoms in Chain A from pdb. Write to a file
    def writeChainACA(self, filename):
        fd=open(filename, 'w')
        count=0
        for i in xrange(0,len(self.atoms)):
            a=self.atoms[i]
            #a.dump()
            if 'CA' in a.name and a.chainID=='A':
                fd.write(a.writeAtom())
                count=count+1
        fd.close()
        if count==0:
            print "No atom written in [%s]!" % (filename)


    # get tip atom list
    #def atomsbytip(self, profile):
    def atomsbytip(self):
        profile = 'AAtips.py'
        cgs=[]
        # loading tip atoms records
        fp=open(profile, 'r')
        lines = fp.readlines()
        fp.close()
        AAtipDict={}
        for i in xrange(0,len(lines)):  
            line = lines[i].strip()
            AAstr = line.split(',')
            AAname = AAstr[0]
            AAtips = AAstr[2]
           
            AAtipDict[AAname]=AAtips.split(' ')
#            print '%s %s' % (AAname, AAtipDict[AAname])

        # iterate all atoms
        #fd=open(filename, 'w') 
        lastAtom = self.atoms[0]
        currentResi = lastAtom.resSeq    
        matchCount=0  
        outputCount=0
        X=0.0
        Y=0.0
        Z=0.0 
        isDone = 0
        for i in xrange(0,len(self.atoms)):
            a = copy.copy(self.atoms[i])
            #if a.chainID!='A':
            #    continue
            if a.resName not in AAtipDict:
                print 'err:%s:: non AA atom %s %d [%s]]' % (self.pdb, a.resName, a.resSeq, a.name.strip())
                continue

            if a.resSeq == currentResi: # if in the same residue
                if isDone == 0:
                    # check against all the tip atoms
                    for tipAtom in AAtipDict[a.resName]:
                        
                        # summing up coordinates if matching
                        if a.name.strip()==tipAtom:
                            
                            X+=a.x
                            Y+=a.y
                            Z+=a.z
                            matchCount+=1
#                            print '%s:: matching %s %d [%s] with [%s], matchCount: %d' % (self.pdb, a.resName, a.resSeq, a.name.strip(), tipAtom, matchCount)
                    # output final coordinates and change flow control variables
                    if matchCount == len(AAtipDict[a.resName]):
#                        print '%s:: matching complete %s %d' % (self.pdb, a.resName, a.resSeq)
                        a.x = X/matchCount
                        a.y = Y/matchCount
                        a.z = Z/matchCount
#                        print '[%s]\n' % a.writeAtom()
                        #fd.write(a.writeAtom())
                        cgs.append(a)
                        outputCount+=1
                        X=0.0
                        Y=0.0
                        Z=0.0 
                        isDone = 1
#                lastAtom = a
#                print 'matchCount trace : %d' % (matchCount)
            else: # residue number changed
                if matchCount!=len(AAtipDict[lastAtom.resName]):
                    #print "%s:: Tip atom not found for [%d] [%s], use last atom [%s] instead. matchCount: %d, tipCount: %d" % (self.pdb, lastAtom.resSeq, lastAtom.resName, lastAtom.name, matchCount, len(AAtipDict[lastAtom.resName]))
                    #fd.write(lastAtom.writeAtom())
                    X=0.0
                    Y=0.0
                    Z=0.0                     
                if a.name.strip()!='N':
                    print "err:%s:: No leading [N] for RESI [%d] [%s]" % (self.pdb, a.resSeq, a.resName)
                currentResi = a.resSeq
                matchCount=0
                isDone=0
            lastAtom = a # save last atom after all business' done
        
        # for the last residue (there is no residue number change for it)
        if matchCount!=len(AAtipDict[lastAtom.resName]):
            print "err:%s:: Tip atom not found for [%d] [%s]" % (self.pdb, lastAtom.resSeq, lastAtom.resName)
            #fd.write(a.writeAtom())
        if outputCount==0:
            print "err:No atom written from %s" % (self.pdb)            

        return cgs

		
    # tip atom extraction
    # not for chain A only
    # def writeChainATips(self, profile, filename):
    # this will change the original atoms!!!!
    def writeTips(self, profile, filename):
        # loading tip atoms records
        fp=open(profile, 'r')
        lines = fp.readlines()
        fp.close()
        AAtipDict={}
        for i in xrange(0,len(lines)):  
            line = lines[i].strip()
            AAstr = line.split(',')
            AAname = AAstr[0]
            AAtips = AAstr[2]
           
            AAtipDict[AAname]=AAtips.split(' ')
#            print '%s %s' % (AAname, AAtipDict[AAname])

        # iterate all atoms
        fd=open(filename, 'w') 
        lastAtom = self.atoms[0]
        currentResi = lastAtom.resSeq    
        matchCount=0  
        outputCount=0
        X=0.0
        Y=0.0
        Z=0.0 
        isDone = 0
        for i in xrange(0,len(self.atoms)):
            #a = self.atoms[i]
            a = copy.copy(self.atoms[i])
            #if a.chainID!='A':
            #    continue
            if a.resName not in AAtipDict:
                print '%s:: non AA atom %s %d [%s]]' % (self.pdb, a.resName, a.resSeq, a.name.strip())
                continue

            if a.resSeq == currentResi: # if in the same residue
                if isDone == 0:
                    # check against all the tip atoms
                    for tipAtom in AAtipDict[a.resName]:
                        
                        # summing up coordinates if matching
                        if a.name.strip()==tipAtom:
                            
                            X+=a.x
                            Y+=a.y
                            Z+=a.z
                            matchCount+=1
#                            print '%s:: matching %s %d [%s] with [%s], matchCount: %d' % (self.pdb, a.resName, a.resSeq, a.name.strip(), tipAtom, matchCount)
                    # output final coordinates and change flow control variables
                    if matchCount == len(AAtipDict[a.resName]):
#                        print '%s:: matching complete %s %d' % (self.pdb, a.resName, a.resSeq)
                        a.x = X/matchCount
                        a.y = Y/matchCount
                        a.z = Z/matchCount
#                        print '[%s]\n' % a.writeAtom()
                        fd.write(a.writeAtom())
                        outputCount+=1
                        X=0.0
                        Y=0.0
                        Z=0.0 
                        isDone = 1
#                lastAtom = a
#                print 'matchCount trace : %d' % (matchCount)
            else: # residue number changed
                if matchCount!=len(AAtipDict[lastAtom.resName]):
                    #print "%s:: Tip atom not found for [%d] [%s], use last atom [%s] instead. matchCount: %d, tipCount: %d" % (self.pdb, lastAtom.resSeq, lastAtom.resName, lastAtom.name, matchCount, len(AAtipDict[lastAtom.resName]))
                    #fd.write(lastAtom.writeAtom())
                    X=0.0
                    Y=0.0
                    Z=0.0                     
                if a.name.strip()!='N':
                    print "%s:: No leading [N] for RESI [%d] [%s]" % (self.pdb, a.resSeq, a.resName)
                currentResi = a.resSeq
                matchCount=0
                isDone=0
            lastAtom = a # save last atom after all business' done
        
        # for the last residue (there is no residue number change for it)
        if matchCount!=len(AAtipDict[lastAtom.resName]):
            print "%s:: Tip atom not found for [%d] [%s]" % (self.pdb, lastAtom.resSeq, lastAtom.resName)
            #fd.write(a.writeAtom())
        if outputCount==0:
            print "No atom written from [%s]!" % (filename)            



    # write concise version of a tip file 
    # x,y,z,resi,resn
    # input tip file 
    # output the concised version
    def writeSPDB(self):
        aamap = AAmap()
        fo = open(self.pdbfile+'.spdb', 'w')
        for a in self.atoms:
            fo.write('%f %f %f %d %s\n' % (a.x, a.y, a.z, a.resSeq, aamap.getAAmap(a.resName)))
        fo.close()
    
    # write all ATOM entries as a pdbfile
    def writepdb(self, filename):
        with open(filename, 'w') as fout:
            fout.write('%s\n' % ''.join([at.writeAtom() for at in atomlist]))

    # use spectral clustering method to find residue contact groups        
    def spectralClustering(self, cluster_size):
        from sklearn.cluster import spectral_clustering
        from scipy.sparse import coo_matrix

        rowlist = []
        collist = []
        datalist = []
        N = len(self.atoms)
        for i in xrange(0,N):
            v1=np.array((self.atoms[i].x, self.atoms[i].y, self.atoms[i].z))
            for j in xrange(0,N):
                v2 = np.array((self.atoms[j].x, self.atoms[j].y, self.atoms[j].z))
                if i == j:
                    euclidean = 0.0
                    affinity = 5.0 # set a large value
                else:
                    euclidean = np.linalg.norm(v1-v2)
                    affinity = 1/euclidean

                key = "%d-%d" % (i,j)
                self.pairwiseDict[key]= euclidean

                rowlist.append(i)
                collist.append(j)
                datalist.append(affinity)

        # prepare affinity matrix for clustering       
        row = np.array(rowlist)
        col = np.array(collist)
        data = np.array(datalist)

        graph = coo_matrix((data, (row, col)), shape=(N, N))
        #labels = spectral_clustering(graph, n_clusters=int(N/cluster_size), eigen_solver='arpack') 
        labels = spectral_clustering(graph, n_clusters=int(N/cluster_size))

        amap = AAmap()
        cluster2fid = {}
        for index, lab in enumerate(labels) :
            cluster2fid.setdefault(lab, [])
            cluster2fid[lab].append(index)

        for key in cluster2fid: # for each cluster
            c=cluster(self.pdb, self.top, self.pfam, '', '', self.seqheader, '', '', self.center, self.cutoff, self.scutoff, self.flag, 1.0, self.desc)
            for index in cluster2fid[key]:
                c.addNeighbor(amap, self.atoms[index],index) 
            c.pdbidx=c.pdbidx.lstrip() # will change meanDist
            c.pdbResSeq=c.pdbResSeq.lstrip()
            meanDist = self.clusterMeanDist(c)
            print ('%s,%0.2f,%s,%s,%s') % (self.pdb, meanDist, ''.join(sorted(c.str)), ''.join(sorted(c.typeStr)), c.pdbResSeq)

    # pairwise
    # read all atom XXXX_A.pdb and get pairwise contact within a cutoff and a sequential distance
    # output XXXX_A.res.csu_d
    def pairContactbyCutoff(self, cutoff, seqdist):
        print self.pdbfile 
        cid = self.pdbfile[5] # just for XXXX_A.pdb naming format
        existDict = {}
        fo = open(self.pdbfile[:-4]+'.res.csu_d', 'w')
        for i in xrange(0, len(self.atoms)):
            for j in xrange(0, len(self.atoms)):
                a = self.atoms[i]
                b = self.atoms[j]
                if abs(a.resSeq - b.resSeq) <= seqdist or a.resSeq == b.resSeq:
                    continue
                v1 = np.array((a.x, a.y, a.z))
                v2 = np.array((b.x, b.y, b.z))
                dist =  np.linalg.norm(v1-v2) 
                key = '%d %d' % (a.resSeq, b.resSeq)
                if key not in existDict and dist <= cutoff:
                    fo.write('%s\t%d%s\t%s\t%d%s\n' % (a.resName, a.resSeq, cid, b.resName, b.resSeq, cid))
                    existDict['%d %d' % (a.resSeq, b.resSeq)] = 1
                    existDict['%d %d' % (b.resSeq, a.resSeq)] = 1
        fo.close()



    # calculate pairwise distance
    def pairwise(self):
        for i in xrange(0,len(self.atoms)):
            v1=np.array((self.atoms[i].x, self.atoms[i].y, self.atoms[i].z))
            for j in xrange(0,len(self.atoms)):
                key = "%d-%d" % (i,j)
                v2 = np.array((self.atoms[j].x, self.atoms[j].y, self.atoms[j].z))
                self.pairwiseDict[key]= np.linalg.norm(v1-v2)
    
    # print pairwise distance between atom[i] and /other atoms   
    # index starts from 0         
    def getPairwiseOf(self, index):
            a = self.atoms[index]
            a.dump()
            v1=np.array((a.x, a.y, a.z))
            for j in xrange(0,len(self.atoms)):
                key = "%d-%d" % (index,j)
                v2 = np.array((self.atoms[j].x, self.atoms[j].y, self.atoms[j].z))
                print "%s: [%f], [%d]" % (key, np.linalg.norm(v1-v2), abs(index-j))
        
    # find clusters (centered by CA). Save all clusters in an array
    def filterClusters(self):
        if len(self.pairwiseDict)==0:
            self.pairwise()
        amap = AAmap()

        for i in xrange(0,len(self.atoms)):
            c=cluster(self.pdb, self.top, self.pfam, '', '', self.seqheader, '', '', self.center, self.cutoff, self.scutoff, self.flag, 1.0, self.desc)
            c.addNeighbor(amap, self.atoms[i],i) # put itself in first
            nbnum=0
            for j in xrange(0,len(self.atoms)):
                key= "%d-%d" % (i, j)
                if (self.pairwiseDict[key] <= self.cutoff) and (abs(i-j) >= self.scutoff):
                    c.addNeighbor(amap, self.atoms[j], j)
                    nbnum=nbnum+1
                    c.thetaPhi.append(self.calculateThetaPhi(self.atoms[i], self.atoms[j]))
            if nbnum<self.nbcutoff: 
                continue
                
            c.pdbidx=c.pdbidx.lstrip() # will change meanDist
            c.pdbResSeq=c.pdbResSeq.lstrip()
            meanDist = self.clusterMeanDist(c)
            if meanDist < 5.8:
                print ('%s,%0.2f,%s,%s,%s,%s') % (self.pdb, meanDist, ''.join(sorted(c.str)), ''.join(sorted(c.typeStr)), c.pdbResSeq, self.getSphericalStr(c))
                self.clusters.append(c)
                #self.writeClusterPDB(('%s%d.pdb') % (self.pdb,len(self.clusters)), templateAtom, c.thetaPhi)                


    def calculateThetaPhi(self, at1, at2):
#        v1 = np.array((at1.x, at1.y, at1.z))
#        v2 = np.array((at2.x, at2.y, at2.z))
#        v0 = (v1-v2)/np.linalg.norm(v1-v2)       
#        
#        x = v0[0]
#        y = v0[1]
#        z = v0[2]
        x=at1.x-at2.x
        y=at1.y-at2.y
        z=at1.z-at2.z
        
        phi = math.atan2(z, math.sqrt(x*x+y*y))
        th = math.atan2(y,x)        
        return (th, phi)
    
    def getSphericalStr(self, c):
        x=''
        y=''
        z=''
        for item in c.thetaPhi:
            x=('%s %0.4f') % (x,item[0])
            y=('%s %0.4f') % (y,item[1])
            z=('%s 0') % (z)
        return (('%s%s%s,%d') % (x,y,z,len(c.thetaPhi))).lstrip()

    # get mean pairwised distance for each atom in the cluster
    # should filter those who has large mean dist value
    # which means they are not real clusters
    def clusterMeanDist(self,cl):
        dist=0.0
        count=0
        idxArray=cl.pdbidx.split(' ')
        for i in xrange(0, len(idxArray)):
            for j in xrange(i+1, len(idxArray)):
                count+=1
                key='%s-%s' % (idxArray[i],idxArray[j])
                dist+=self.pairwiseDict[key]
        #print 'clusterMeanDist:: ', dist/count
        if count==0:
            return 0.0
        return dist/count

                    
    # print a cluster object
    def dumpClusters(self):
        for i in xrange(0,len(self.clusters)):
            a=self.atoms[i]
            c=self.clusters[i]
            print '++++++++++++++++++++++++++++++++++\n'
            a.dump()
            c.dump()
            print '++++++++++++++++++++++++++++++++++\n' 
    
    # print protein object
    def dump(self):
 		print ('pdb:[%s]\n' + 
		 	  'top:[%s]\n' +
		 	  'pfam:[%s]\n' +
		 	  'center:[%s]\n'+      
		 	  'cutoff:[%f]\n' +
		 	  'scutoff:[%d]\n' +
		 	  'seqheader:[%s]\n' +
		 	  'flag:[%d]\n' +
		 	  'desc:[%s]\n') % \
		 	  (
		 	  	self.pdb,
		 	  	self.top,
		 	  	self.pfam,
		 	  	self.center,
		 	  	self.cutoff,
		 	  	self.scutoff,
		 	  	self.seqheader,
		 	  	self.flag,
		 	  	self.desc
		 	  )

    def getClusterNum(self):
		return len(self.clusters) 

    # return sliced pdb atom array
    # s: input seq(substring) to match pdb seq
    def slicebyseq(self, s):
        pdbseq = self.seq.upper()
        idx = pdbseq.find(s)
        return self.resAtoms[idx:idx+len(s)] if idx != -1 else []



