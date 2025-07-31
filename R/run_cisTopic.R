source('cisTopic.R')
cores = n_seeds_cisTopic
cl <- makeCluster(cores) #not to overload your computer
registerDoParallel(cl)
ari_scores = rep(NA,n_seeds_cisTopic)

ari_scores <- foreach(i=1:n_seeds_cisTopic,.combine=append) %dopar% {
  suppressWarnings(library(cisTopic))
  require(aricode)
  pred = cisTopic_fn(cisTopicObject,i,n_ct,res_pth)
  ari_score = ARI(labels[names(pred)],pred)
  ari_score
  #do other things if you want
}
score = median(ari_scores)
write(score,file=paste(res_pth,'ari_score.txt',sep=''))
stopCluster(cl)


