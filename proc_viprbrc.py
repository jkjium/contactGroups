import commp as cp

def splitbyprotein(arglist):
    if len(arglist) < 1:
        cp._err('Usage: python proc_viprbrc.py splitbyprotein 93581267956-ProteinFastaResults.fasta')
    
    infile = arglist[0]
    namedict = {
               'envelope_protein': 'e_protein',
               'e_protein':'e_protein',
               'rna_dependent_rna_polymerase':'rna_dependent_rna_polymerase',
               'leader_protein': 'leader_protein',
               'nucleocapsid_phosphoprotein':'n_protein',
               'nucleocapsid_protein':'n_protein',
               'n_protein': 'n_protein',
               'nucleoprotein': 'n_protein',
               'spike_glycoprotein':'s_protein',
               'spike_protein':'s_protein',
               'surface_glycoprotein': 's_protein',
               's_protein': 's_protein',
               '3c_like_proteinase': '3c_like_proteinase',
               'endornase': 'endornase',
               '3\'_to_5\'_exonuclease':'3_to_5_exonuclease',
               'membrane_glycoprotein':'m_protein',
               'membrane_protein': 'm_protein',
               'm_protein': 'm_protein',
               'helicase': 'helicase'
    }

    # open file handles
    filedict = {}
    for f in set(namedict.values()):
        fp = open(f+'.fa', 'w')
        filedict[f] = fp
        
    # split sequence into different files
    for h, s in cp.fasta_iter(infile):
        sarr = h.split('|')
        proteintype = sarr[1].lower()
        if proteintype in namedict:
            filedict[namedict[proteintype]].write('>%s\n%s\n' % (h, s))
        
    # close file
    for k in filedict:
        filedict[k].close()
    
    cp._info('sequence split done.')

if __name__=='__main__':
    cp.dispatch(__name__)