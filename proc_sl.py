import commp_ns as cp
import numpy as np
import json

# extract eggnog id by input
# eggstr: Eukaryota,38C7D@33154
# group_name: Eukaryota
def _eid(eggstr, group_name):
    ret = 'na na'
    for t in eggstr.split('|'):
        if group_name in t:
            t1 = t.split(',')
            ret = '%s %s' % (t1[0], t1[1]) if len(t1) == 2 else '%s null' % (t)
            break
    return ret


# input: kjia@DESKTOP-DKBVSMN /cygdrive/d/workspace/sponge/embedding/sl/annotations/eggnog
# extract eggnog information from eggnog.4.sl.emb.tsv
# c73883_g1_i1_m.5363	AF-Q2F637-F1	COG5040@1|root,KOG0841@2759|Eukaryota,38CKG@33154|Opisthokonta,3BCFK@33208|Metazoa,3CT6T@33213|Bilateria,41U1K@6656|Arthropoda,3SKBS@50557|Insecta,46KNG@7399|Hymenoptera	c73883_g1_i1_m.5363	Q04917	sp|Q04917|1433F_HUMAN	COG5040@1|root,KOG0841@2759|Eukaryota,38CKG@33154|Opisthokonta,3BDYD@33208|Metazoa,3D3HT@33213|Bilateria,48A1P@7711|Chordata,4955S@7742|Vertebrata,3J48D@40674|Mammalia,35AKZ@314146|Euarchontoglires,4MH1E@9443|Primates,35Z69@314294|Cercopithecoidea
def eggnog_cmp(args):
    assert len(args) == 2, 'Usage: python proc_sl.py eggnog_cmp eggnog.4.sl.emb.tsv out'
    infile = args[0]
    outfile = args[1]

    eid_common_euk=0
    eid_common_met=0
    eid_common_euk_met=0

    outlist =[]
    for s in cp.loadtuples(infile, '\t'):
        n = s[0] # protein name
        fs = s[2] # foldseek result
        es = s[6] # emb result

        f_eid_e = _eid(fs, 'Eukaryota')
        f_eid_m = _eid(fs, 'Metazoa')

        e_eid_e = _eid(es, 'Eukaryota')
        e_eid_m = _eid(es, 'Metazoa')

        outlist.append('%s %s %s %s %s' % (n,f_eid_e,f_eid_m,e_eid_e,e_eid_m))

        # count # of identical Eukaryota
        if f_eid_e == e_eid_e:
            eid_common_euk+=1
        if f_eid_m == e_eid_m:
            eid_common_met+=1
        if (f_eid_e == e_eid_e) and (f_eid_m == e_eid_m):
            eid_common_euk_met+=1

    print('eid_common_euk %d' % eid_common_euk)
    print('eid_common_met %d' % eid_common_met)
    print('eid_common_euk_met %d' % eid_common_euk_met)

    with open(outfile, 'w') as fout:
        fout.write('%s\n' % '\n'.join(outlist))
    cp._info('save sbs id extraction to %s' % outfile)


# parse json file from uniprot curl queries
# curl "https://rest.uniprot.org/uniprotkb/search?query=A0A024B5K5&fields=lineage"
# input: json list file; each line is a json result
# output: outprefix.uniprot.lineage.txt
# output: outprefix.sname.taxid.rank.flag.txt
def parse_taxjson(args):
    assert len(args) == 2, 'Usage: python proc_sl.py parse_taxjson jsonfile outprefix'
    jsonfile = args[0]
    outprefix = args[1]

    taxdb = set()
    outlist = []
    for line in cp.loadlines(jsonfile):
        j = json.loads(line)
        uniprotid = j['results'][0]['primaryAccession']
        if 'lineages' not in j['results'][0]:
            continue
        lineages = j['results'][0]['lineages']
        lineage_list = []
        for e in lineages:
            scientificname = e['scientificName']
            #print(scientificname)
            #commonname = e['commonName']
            taxid = e['taxonId']
            rank = e['rank']
            taxdb.add('%d,%s,%s' % (taxid, scientificname, rank))
            lineage_list.append(scientificname)
        lineage_list.reverse()
        outlist.append('%s|%s' % (uniprotid, ','.join(lineage_list)))

    # output lineage series 
    outfile = '%s-uniprot-lineage.txt' % outprefix
    with open(outfile, 'w') as fout:
        fout.write('%s\n' % '\n'.join(outlist))
    cp._info('save lineage series to %s' % outfile)

    # output tax db
    outfile = '%s-tax-sname-rank.txt' % outprefix
    with open(outfile, 'w') as fout:
        fout.write('%s\n' % '\n'.join(list(taxdb)))
    cp._info('save tax db to %s' % outfile)


# calculate statistics:
# generate side-by-side tsv for sl proteni annotations
# input: cmp-2-uniprot-lineage.txt
def lineage_cmb(args):
    assert len(args) == 3, 'Usage: python proc_sl.py lineage_cmp cmp-1-sl-em-common.vec2 cmp-2-uniprot-lineage.txt outfile'
    cmpfile = args[0]
    linfile = args[1]
    outfile = args[2]

    lineagedb = dict((v[0], v[1]) for v in cp.loadtuples(linfile,'|'))
    cp._info('%d lineage info loaded.' % len(lineagedb))

    outlist = []
    for c in cp.loadtuples(args[0]):
        sp = c[0]
        fs_uniprot = c[1]
        pr_uniprot = c[2]
        fs_anno = lineagedb[fs_uniprot] if fs_uniprot in lineagedb else 'obsolete'
        pr_anno = lineagedb[pr_uniprot]
        outlist.append('%s\t%s\t%s\t%s\t%s' % (sp, fs_uniprot, fs_anno, pr_uniprot, pr_anno))

    with open(outfile, 'w') as fout:
        fout.write('%s\n' % '\n'.join(outlist))
    cp._info('save to %s' % outfile)


# calculating identical annotation by different branch levels
# input: cmp-3-sbs-sl-em-anno.tsv
# output: number of identicals
def lineage_cmp(args):
    assert len(args) == 2, 'Usage: python proc_sl.py lineage_cmp cmp-3-sbs-sl-em-anno.tsv level'
    infile = args[0]
    n = int(args[1])
    sum = 0
    retstr = set()
    for line in cp.loadtuples(infile,'\t'):
        fs_anno = line[2].split(',')
        pr_anno = line[4].split(',')
        fs_str = ','.join(fs_anno[:n])
        pr_str = ','.join(pr_anno[:n])
        if fs_str==pr_str:
            sum+=1
            retstr.add(pr_str)
    print('%d\n%s\n' % (sum, '\n'.join(retstr)))

# extract the best hit from prost results
def besthit_prost(args):
    assert len(args) == 2, 'Usage: python proc_sl.py besthit_prost run2.h.tsv run2.h.prost.best.tsv'
    infile = args[0]
    outfile = args[1]

    besthits = set() 
    outlist = []
    for s in cp.loadtuples(infile, '\t'):
        # {sname, uniprotID, prost_dist, evalue}
        if s[0] not in besthits:
            outlist.append('%s,%s,%s,%s' % (s[0],s[1],s[5],s[6]))
            besthits.add(s[0])

    with open(outfile, 'w') as fout:
        fout.write('%s\n' % '\n'.join(outlist))
    cp._info('save prost best hits to %s' % outfile)


if __name__=='__main__':
    cp.dispatch(__name__)







