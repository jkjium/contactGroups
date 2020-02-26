import json
import commp as cp

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
        # if within disordered region
        if (msai1 not in m) or (msai2 not in m):
            continue
        if m[msai1] >= dstart and m[msai2] <= dend:
            outstr.append('%s D %s %d %d' % (line, disprotid, m[msai1], m[msai2]))
            dcount+=1
        else:
            outstr.append('%s O %s %d %d' % (line, disprotid, m[msai1], m[msai2]))
            ocount+=1
    
    # output a new .dtype.cflat file
    with open(outfile, 'w') as fout:
        fout.write('%s\n' % '\n'.join(outstr))
    print ('%s %d %d %d %d' % (outfile, plen, dlen, dcount, ocount))




if __name__=='__main__':
    cp.dispatch(__name__)