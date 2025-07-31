#jaspar_motifs = readJASPARMatrix(paste(workspace_dir_pth,'jaspar_motifs/',tolower(org),'.txt',sep=''),matrixClass=c('PFM'))
motif_ix <- matchKmers(6, fragment_counts, 
                      genome = bsgenome)
#rm(jaspar_motifs)
