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


if __name__=='__main__':
    cp.dispatch(__name__)