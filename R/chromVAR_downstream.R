sample_cor <- cor(t(cell_motif), use = "pairwise.complete.obs")
save(sample_cor,file=paste(res_pth,'sample_cor.RData',sep=''))
rm(cell_motif)
clust = pheatmap(as.dist(sample_cor),
                 clustering_distance_rows = as.dist(1-sample_cor),
                 clustering_distance_cols = as.dist(1-sample_cor))
rm(sample_cor)
clust = clust$tree_row
save(clust,file=paste(res_pth,'clust.RData',sep=''))
pred = cutree(clust,n_ct)
rm(clust)
save(pred,file=paste(res_pth,'pred.RData',sep=''))
score = ARI(pred,labels[names(pred)])
write(score,file=paste(res_pth,'ari_score.txt',sep=''))
