import json
import commp as cp
import collections

# from disprot.json extract all the entries that pfam region contains disordered region
# output: disportid pdb pfamid pfam_start:end disordered_start:end %_of_disorder pfam_seq
def filterbypfam(arglist):
    if len(arglist) < 2:
        cp._err('Usage: python proc_disprot.py filterbypfam disprot.json outfile') 
    # t.json
    jsonfile = arglist[0]
    outfile = arglist[1]
    with open(jsonfile) as fp:
        data = fp.read()
    raw = json.loads(data)
    d = raw['data']
    '''
    # try to extract pdb id
    for r in x['regions']:
        if 'cross_refs' in r:
            for c in r['cross_refs']:
                if c['db'] == 'PDB':
                    print c['id']
    '''
    outlist = []
    count = 0
    for x in d:
        count+=1
        disprotid =  x['disprot_id']
        seq = x['sequence']
        pdb = 'na'
        # try to extract pdb id
        for r in x['regions']:
            if 'cross_refs' in r:
                for c in r['cross_refs']:
                    if c['db'] == 'PDB':
                        pdb = c['id']
        pfamlist = []
        disorderlist = []
        if 'pfam' in x['features']:
            pfam = x['features']['pfam']
            # p[0]:start, p[1]:end, p[2]: id
            for p in  pfam:
                pfamlist.append((p['start'],p['end'], p['id']))
            # d[0]:start d[1]:end
            for r in x['disprot_consensus']['structural_state']:
                if r['type'] == 'D':
                    disorderlist.append((r['start'], r['end']))
        
        if (len(pfamlist)!=0) and (len(disorderlist)!=0):
            #print disprotid, len(pfamlist), len(disorderlist)
            for p in pfamlist:
                for d in disorderlist:
                    # just in case
                    pstart = int(p[0])
                    pend = int(p[1])

                    dstart = int(d[0])
                    dend = int(d[1])
                    # filter pfam region contains disordered region
                    if (pstart <= dstart) and (pend >= dend):
                        outlist.append('%s %s %s %d %d %d %d %.2f %s' % (disprotid, pdb, p[2], pstart, pend, dstart, dend, 1.0*(dend-dstart)/(pend-pstart), seq))

    with open(outfile, 'w') as fout:
        fout.write('%s\n' % '\n'.join(outlist))
    cp._info('%d / %d entries are saved to %s' % (len(outlist), count, outfile))



# .dtype.cflat format
# 0.pdbid 1.chainID 2.r1 3.r2 4.rn1 5.rn2 6.sc.d 7.ca.d 8.tip.d 9.area1 10.area2 11.PfamID 12.m1 13.m2 14.dca 15.trans_dca 16.trans_dca_zscore 17.mi 18.trans_mi 19.trans_mi_zscore 20.mip 21.trans_mip 22.trans_mip_zscore	
# 23.dtype 24.disprotid 25.orig_posi1 26.orig_posi2	
# origpos: position id from the original disport sequence
# scan through each amino acid
# sum up ce value with other properties (specified in "opt")
def scancflat(arglist):
    if len(arglist) < 3:
        cp._err('Usage: python proc_disprot.py scancflat infile.cflat opt outfile')
    infile = arglist[0]
    # DP00021-PF00009-10-282-46-67-251-21-dtype.cflat
    sarr0 = infile.split('-')
    dstart = int(sarr0[4])
    dend = int(sarr0[5])


    opt = int(arglist[1])
    outfile = arglist[2]

    outdict = collections.defaultdict(float)
    orderdict = collections.defaultdict(float)
    disdict = collections.defaultdict(float)
    corder=0
    cdis=0
    for line in cp.loadlines(infile):
        sarr = line.split(' ')
        pos1 = int(sarr[25])
        pos2 = int(sarr[26])
        if sarr[14] == '-191':
            continue
        ce = float(sarr[14]) # dca
        score = 0 
        if opt == 0:
            score = ce
        elif opt == 1:
            score = ce * (pos2 - pos1)
        elif opt == 2:
            if float(sarr[7]) < 5 and sarr[16] > 0.0:
                score =  ce
        elif opt == 3:
            score = ce * (pos2 - pos1) * (pos2 - pos1)
        elif opt == 4:
            if float(sarr[7]) < 5:
                score = ce

        outdict[pos1]+=score
        outdict[pos2]+=score

        # interaction within disordered region only
        if (pos1>=dstart) and (pos2<=dend):
            disdict[pos1]+=score
            disdict[pos2]+=score
            cdis+=1
        if (pos1>dend) or (pos2<dstart): # within the ordered region before / after the disordered 
            orderdict[pos1]+=score
            orderdict[pos2]+=score
            corder+=1
        if (pos1<dstart) and (pos2>dend): # interaction between ordered regions before and after the disordered
            orderdict[pos1]+=score
            orderdict[pos2]+=score
            corder+=1


        '''
        # interactions within disorder and interactions from disordered region with ordered region
        if (pos1<dstart or pos1>dend) and (pos2<dstart or pos2>dend):
            orderdict[pos1]+=score
            orderdict[pos2]+=score
            corder+=1
        else: # interactions ordered region only
            if (pos1>=dstart) and (pos1<=dend): 
                disdict[pos1]+=score
            if (pos2>=dstart) and (pos2<=dend):
                disdict[pos2]+=score
            cdis+=1
        '''

    #print '%s %.4f %4f %d %d %d %d' % (infile, sum(orderdict.values())/corder, sum(disdict.values())/cdis, len(orderdict), len(disdict), corder,cdis)
    print '%s %.4f %4f %d %d %d %d' % (infile, sum(orderdict.values())/corder, sum(disdict.values())/cdis, len(orderdict), len(disdict), dstart,dend)

    total = sum(outdict.values())
    x = outdict.keys()
    x.sort()
    with open(outfile, 'w') as fout:
        fout.write('%s\n' % ('\n'.join(['%d %.4f' % (i, outdict[i]/total) for i in x])))
        #fout.write('%s\n' % ('\n'.join(['%d %.4f' % (i, outdict[i]) for i in x])))


# extract cflat2 records by disordered / non-disordered regions
# output two files for each disordered record: .disorder.cflat, .order.cflat
def splitcflat2(arglist):
    if len(arglist) < 2:
        cp._err('Usage: python proc_disprot.py splitcflat2 seq2msamap.map PF00000.cflat2')
    
    mapfile = arglist[0]
    cflatfile = arglist[1]

    # extracting pfam, region indices info from mapfile name
    # DP00822-PF00472-7-132-100-130.fa-PF00472.map
    sarr = mapfile.split('.')
    # DP00822-PF00472-7-132-100-130
    disprotid = sarr[0]
    # output files

    sarr = disprotid.split('-')
    pfam = sarr[1]
    # pfam region
    pstart = int(sarr[2])
    pend = int(sarr[3])
    # disordered region
    dstart = int(sarr[4])
    dend = int(sarr[5])

    dlen = dend - dstart
    plen = pend - pstart - dlen

    outfile = '%s-%d-%d-dtype.cflat' % (disprotid, plen, dlen)

    # a dictionary from msai to posi 
    # map index
    # >DP00822-PF00472-7-132-100-130
    # NVHLPDAEIELTAIRAQGAGGQNVNKVSSAMHLRFDINASSLPPFYKERLLALNDSRITSDGVIVLKAQQYRTQEQNRADALLRLSELIVNAAKVEKKRRPTRPTLGSKTRRLESKSKRGSIKAGR
    # N:0, 7 (not in map)
    # V:1, 8
    # H:2, 9
    # map resi +7 = the original index
    # in map file:
    # posi posn msai msan
    # 1 V 75 V
    # 2 H 76 H
    m = {}
    for line in cp.loadlines(mapfile):
        sarr = line.split(' ')
        m[int(sarr[2])] = int(sarr[0])+pstart # add pstart to restore the original index from disprot db

    # .cflat2 format
    # 0.pdbid 1.chainID 2.r1 3.r2 4.rn1 5.rn2 6.sc.d 7.ca.d 8.tip.d 9.area1 10.area2 11.PfamID 12.m1 13.m2 
    # 14.dca 15.trans_dca 16.trans_dca_zscore 17.mi 18.trans_mi 19.trans_mi_zscore 20.mip 21.trans_mip 22.trans_mip_zscore

    # append type 'O' or 'D', disprotid, orig_pos1, orig_pos2 to cflat2
    outstr = []
    dcount = 0
    ocount = 0
    for line in cp.loadlines(cflatfile):
        sarr = line.split(' ')
        msai1 = int(sarr[12])
        msai2 = int(sarr[13])
        if (msai1 not in m) or (msai2 not in m):
            continue
        # if pair within disordered region, mark pair as type "D"
        if m[msai1] >= dstart and m[msai2] <= dend:
            outstr.append('%s D %s %d %d' % (line, disprotid, m[msai1], m[msai2]))
            dcount+=1
        else: # mark as "O"
            outstr.append('%s O %s %d %d' % (line, disprotid, m[msai1], m[msai2]))
            ocount+=1
    
    # output a new .dtype.cflat file
    with open(outfile, 'w') as fout:
        fout.write('%s\n' % '\n'.join(outstr))
    print ('%s %d %d %d %d' % (outfile, plen, dlen, dcount, ocount))


# window slide smoothing to detect disordered regions in the sequence
# call after scancflat
# infile .cflat.scan
# output: .scan.ws
def windowslide(arglist):
    if len(arglist) < 3:
        cp._err('Usage: python proc_disprot.py windowslide infile.scan window_size outfile.ws')

    infile = arglist[0]
    wsize = int(arglist[1])
    outfile = arglist[2]

    scorelist = []
    ticklist = []
    for line in cp.loadlines(infile):
        sarr = line.split(' ')
        ticklist.append(sarr[0])
        scorelist.append(float(sarr[1]))
    
    #smlist = [sum(scorelist[i:wsize]) for i in xrange(0, scorelist)]
    smlist = []
    fout = open(outfile, 'w')
    for i in xrange(0, len(scorelist)):
        slide = scorelist[i:i+wsize]
        smlist.append(sum(slide)/len(slide))
        fout.write('%s %.8f\n' % (ticklist[i], smlist[i]))
    fout.close()
    cp._info('save to %s' % outfile)

if __name__=='__main__':
    cp.dispatch(__name__)