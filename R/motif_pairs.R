load(paste(chromVAR_org_pth,'motif_presence_FIMO_counts_dgC.RData',sep=''))
m_p = as(annotationMatches(motif_ix),'sparseMatrix')
motifs = colnames(m_p)

n = dim(m_p)[2]
res = list()
count = 1
if(same){
 for (i in 1:n){
  for (j in i:n){
      res[[count]] = c(i,j)
      count = (count+1)
  }
}
}else{
 for (i in 1:(n-1)){
  for (j in (i+1):n){
      res[[count]] = c(i,j)
      count = (count+1)
  }
}
}  

n_pairs = length(res)

test = list()
motif_pairs = rep(NA,n_pairs)

for (i in 1:n_pairs){
  vec = res[[i]]
  if(vec[1]==vec[2]){
    v = m_p_counts[,vec[1]]
    v[which(v>1)] <- 1
    test[[i]] = v
  }
  else{
    test[[i]] = m_p[,vec[1]]&m_p[,vec[2]]+0
  }
  motif_pairs[i] = paste(motifs[vec[1]],motifs[vec[2]],sep='.')
  rm(vec)
}
rm(res)
test <- lapply(test, as, "sparseMatrix")
m_p_pairs = do.call(cbind,test)
rm(test)
rownames(m_p_pairs) = rownames(m_p)
colnames(m_p_pairs) = motif_pairs
rm(motif_pairs)

m_p = as(cbind(m_p,m_p_pairs),'sparseMatrix')
rm(m_p_pairs)

source('convert_motif_ix_chromVAR.R')
rm(m_p)
