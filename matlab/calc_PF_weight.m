
msafile = '1k2p_PF07714_full.fa'
W = evfold_weight(msafile, '1k2p', 0.3);
size(W)
dlmwrite(strcat(msafile, '.weight'), W, 'delimiter', ',', 'precision', '%10.8f');

