cisTopic_fn <- function(cisTopicObject=NULL,seed_val=NULL,n_ct=NULL,res_pth=NULL){
  stopifnot(inherits(cisTopicObject,'cisTopic'))
  stopifnot(inherits(seed_val,'numeric') | inherits(seed_val,'integer'))
  stopifnot(inherits(res_pth,'character'))
  stopifnot(inherits(n_ct,'numeric') | inherits(n_ct,'integer'))
  
  seed = 987+(seed_val-1)
  cisTopicObject <- runCGSModels(cisTopicObject, topic=50, seed=seed, nCores=1, burnin = 500, iterations = 600, addModels=FALSE)
  cisTopicObject <- selectModel(cisTopicObject,type='maximum')
  cisTopicObject <- runtSNE(cisTopicObject, target='cell',check_duplicates=FALSE)
  save(cisTopicObject,file=paste(res_pth,'cisTopicObject_',seed,'.RData',sep=''))
  clust_sne = cisTopicObject@dr$cell$tSNE
  rm(cisTopicObject)
  ward <- hclust(dist(clust_sne, method='euclidean'),method = 'ward.D')
  rm(clust_sne)
  pred = cutree(ward,n_ct)
  rm(ward)
  pred_pth = paste(res_pth,'pred/',sep='')
  dir.create(pred_pth,showWarnings = FALSE)
  #system(paste('mkdir',pred_pth))
  save(pred,file=paste(pred_pth,'pred_',seed,'.RData',sep=''))
  return (pred)
}

