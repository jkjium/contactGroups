import commp as cp
import numpy as np

def tt(arglist):
    # for test rountine
    print 'hello'
    pass

# read MSA.fa
# input target header (human P53)
# output the number of alphabet(amino acid) difference with the input sequence
def hamming_diff_with(arglist):
    if len(arglist) <2:
        cp._err('Usage: python proc_pf00870.pf hamming_diff_with msafile.fa target_header')
    msafile = arglist[0]
    target = arglist[1]
    targetseq = ''

    msalist = []
    for head, seq in cp.fasta_iter(msafile):
        if target in head:
            targetseq = seq
        msalist.append((head, seq))
    if targetseq == '':
        cp._err('target %s not found' % target)
    for head, seq in msalist:
        diff = [(str(i),targetseq[i],seq[i]) for i in xrange(0, len(targetseq)) if targetseq[i]!=seq[i]]
        diff_idx = ','.join(d[0] for d in diff)
        target_seg = ''.join([d[1] for d in diff])
        seq_seg = ''.join(d[2] for d in diff)
        print '%s %d %s %s %s' % (head, len(diff), target_seg, seq_seg, diff_idx)

if __name__=='__main__':
    cp.dispatch(__name__)