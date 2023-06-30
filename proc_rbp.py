import commp_ as cp
import json

# find field key words here: https://www.uniprot.org/help/return_fields
# curl "https://rest.uniprot.org/uniprotkb/search?query=O60870&fields=gene_names,ft_domain,ft_repeat,ft_motif,ft_region,ft_zn_fing"
# input: p-info.json.list, each line is a json record containing protein information
# input: motif.txt, mapping between original data motif name to uniprot motif name
# output: uniprot_ID,gene_name,RNA_binding_motif,motif_type,keyword_type,start,end
def parse_json(args):
    assert len(args) == 3, 'Usage: python proc_rbp.py parse_json p-info.json.list motif.txt outfile'
    jsonfile = args[0]
    mfile = args[1]
    outfile = args[2]

    # check whether RNA binding motif keywords occurring in motif descriptions
    # RRM is in RRM 1 
    def isMotif(mstr,motifs):
        for m in motifs:
            if m in mstr:
                return m
        return False

    # load motifs to scan
    motifs = set([line[1] for line in cp.loadtuples(mfile, delimiter=',')])

    # scan for rna binding motifs
    outlist = []
    for line in cp.loadlines(jsonfile):
        j = json.loads(line)
        uniprot_id = j['results'][0]['primaryAccession']
        gene_name  = j['results'][0]['genes'][0]['geneName']['value']

        for f in j['results'][0]['features']:
            mtype = f['type']
            motif = f['description']
            #print(mtype, motif)
            # zinc finger is a type, rest motifs are descriptions
            ktype = isMotif(motif, motifs)
            if mtype == "Zinc finger":
                ktype = "Zinc"
                mstart = f['location']['start']['value']
                mend = f['location']['end']['value']
                outlist.append('%s,%s,%s,%s,%s,%s,%s' % (uniprot_id,gene_name,ktype,mtype,motif,mstart,mend))
            elif ktype!=False:
                mstart = f['location']['start']['value']
                mend = f['location']['end']['value']
                outlist.append('%s,%s,%s,%s,%s,%s,%s' % (uniprot_id,gene_name,ktype,mtype,motif,mstart,mend))

    # each line for one motif
    with open(outfile, 'w') as fout:
        fout.write('%s\n' % '\n'.join(outlist))
    cp._info('save to %s' % outfile)


def foo(args):
    cp._info(args)
    pass

if __name__=='__main__':
    cp.dispatch(__name__)
