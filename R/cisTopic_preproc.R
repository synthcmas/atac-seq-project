mat <- cell_motif
rm(cell_motif)
col.med <- apply(X = mat, MARGIN = 2, FUN = median)

dat = t(replicate(dim(mat)[1], col.med))
t_idx = which(mat >= dat)
f_idx = which(mat < dat)
mat[t_idx] = 1
mat[f_idx] = 0

fake_cols = c()
for (i in 1:dim(mat)[2]){
  s = i*1000
  fake_cols = append(fake_cols,paste('chrA:',s,'-',s+500,sep=''))
}
colnames(mat) = fake_cols

print(class(mat))
mat <- as(mat, "sparseMatrix") 
mat = t(mat)
