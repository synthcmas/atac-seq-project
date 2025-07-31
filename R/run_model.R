workspace_dir_pth = '/home/gaurik/nnp_proj/'
data_dir = 'data/furlong/'
data_dir_pth = paste(workspace_dir_pth,data_dir,sep='')
R_dir_pth = paste(workspace_dir_pth,'R/',sep='')
setwd(R_dir_pth)

cisTopic=FALSE
same = FALSE
n_seeds_cisTopic=8
binary=FALSE
chromVAR_preproc=FALSE
FIMO=TRUE
n_cores = 40
withPairs=FALSE
chromVAR_norm=TRUE
bias_correct=TRUE
chromVAR_postproc=TRUE
chromVAR_downstream=FALSE
topic_modeling=FALSE

stopifnot(inherits(data_dir_pth,'character'))
stopifnot(inherits(cisTopic,'logical'))
stopifnot(inherits(binary,'logical'))
stopifnot(inherits(chromVAR_preproc,'logical'))
stopifnot(inherits(FIMO,'logical'))
stopifnot(inherits(withPairs,'logical'))
stopifnot(inherits(chromVAR_norm,'logical'))
stopifnot(inherits(bias_correct,'logical'))
stopifnot(inherits(chromVAR_downstream,'logical'))
stopifnot(inherits(topic_modeling,'logical'))
stopifnot(inherits(chromVAR_postproc,'logical'))
stopifnot(!(chromVAR_downstream & topic_modeling))
if (cisTopic | topic_modeling){
  stopifnot((inherits(n_seeds_cisTopic,'numeric') | inherits(n_seeds_cisTopic,'integer'))  & (n_seeds_cisTopic >= 5))
  n_cores = n_seeds_cisTopic
}else{
  n_cores = n_cores
}

dataset_pth = paste(data_dir_pth,'dataset.RData',sep='')
name_dir = tail(strsplit(data_dir_pth,split='/')[[1]],n=1)
labels_pth = paste(data_dir_pth,'labels.RData',sep='')
bed_file_pth = paste(data_dir_pth,'bed_file.bed',sep='')
con=file(paste(data_dir_pth,'genome_org.txt',sep=''),open="r")
line=readLines(con) 
bsgenome = line[1]
bsgenome_split = strsplit(bsgenome,split='.',fixed = TRUE)[[1]]
genome = bsgenome_split[4]
org = line[2]
close(con)
load(labels_pth)
n_ct = length(unique(labels))

results_main_dir = paste(workspace_dir_pth,'results/',sep='')
dir.create(results_main_dir,showWarnings = FALSE)

results_dir = paste(results_main_dir,name_dir,'/',sep='') 
dir.create(results_dir,showWarnings = FALSE)

proc = paste('cisTopic',cisTopic,'n_seeds_cisTopic',n_seeds_cisTopic,
               'binary',binary,'chromVAR_preproc',chromVAR_preproc,'FIMO',FIMO,
               'same',same,'withPairs',withPairs,'chromVAR_norm',chromVAR_norm,
               'bias_correct',bias_correct,'chromVAR_postproc',
               chromVAR_postproc,'chromVAR_downstream',chromVAR_downstream,
               'topic_modeling',topic_modeling,sep='_')
res_pth = paste(results_dir,proc,'/',sep='')
dir.create(res_pth,showWarnings = FALSE)

chromVAR_pth = paste(res_pth,'chromVAR/',sep='')
dir.create(chromVAR_pth,showWarnings = FALSE)

chromVAR_org_pth = paste(chromVAR_pth,tolower(org),'/',sep='') 
dir.create(chromVAR_org_pth,showWarnings = FALSE)

source('load_libraries.R')

if(cisTopic){
  load(dataset_pth)
  stopifnot(inherits(mat,'dgCMatrix'))
  
  res_pth = paste(results_dir,'cisTopic','/',sep='')
  dir.create(res_pth,showWarnings = FALSE)
  
  source('create_cisTopicObject.R')
  source('run_cisTopic.R')
}else{
  prmpt = paste('fragment_counts_binary_',binary,'.RData',sep='')

  load(dataset_pth)
  stopifnot(inherits(mat,'dgCMatrix'))
  if(binary){
    mat[which(mat!=0)] <- 1
  }
  mat[,which(colSums(mat)!=0)]
  source('chromVAR_setup.R')
  save(fragment_counts,file=paste(chromVAR_pth,prmpt,sep=''))
  
  print('ChromVAR Setup Step Done!')
  
  prmpt = paste('fragment_counts_binary_',binary,'_addBias_',bias_correct,'.RData',sep='')
  
  if(bias_correct){
    source('chromVAR_addBias.R')
  }
  save(fragment_counts,file=paste(chromVAR_pth,prmpt,sep=''))
  
  print('Bias Step Done!')
  
  if ('lib_size.RData'%in%list.files(data_dir_pth)){
    lib_size_pth = paste(data_dir_pth,'lib_size.RData',sep='')
  }else{
    chromVAR_preproc=FALSE
  }
  
  prmpt = paste('fragment_counts_binary_',binary,'_addBias_',bias_correct,'_chromVAR_preproc_',chromVAR_preproc,'.RData',sep='')
  
  if(chromVAR_preproc){
    source('chromVAR_preprocess.R')
  }
  barcodes = data.frame(barcodes=colnames(fragment_counts))
  write.csv(barcodes,file=paste(res_pth,'barcodes.csv',sep=''),quote=FALSE,row.names=FALSE)
  l = data.frame(labels=labels[barcodes$barcodes])
  write.csv(l,file=paste(res_pth,'labels.csv',sep=''),quote=FALSE,row.names=FALSE)
  save(fragment_counts,file=paste(chromVAR_pth,prmpt,sep=''))
  
  print('Preprocessing Step done!')
  
  prmpt = paste('motif_ix_FIMO_',FIMO,'.RData',sep='')
  
  if(FIMO){
    source('find_motif_FIMO.R')
  }else{
    source('find_motif_chromVAR.R')
  }
  save(motif_ix,file=paste(chromVAR_org_pth,prmpt,sep=''))
  
  print('Motif finding Step Done!')
  
  prmpt = paste('motif_ix_FIMO_',FIMO,'_same_',same,'_withPairs_',withPairs,'.RData',sep='')
  
  if(withPairs){
    source('motif_pairs.R')
  }
  save(motif_ix,file=paste(chromVAR_org_pth,prmpt,sep=''))
  
  print('Motif pairs Step Done!')
  
  prmpt_dev = paste('dev_chromVAR_norm_',chromVAR_norm,'_FIMO_',FIMO,'_same_',same,'_withPairs_',withPairs,'_binary_',binary,'_addBias_',bias_correct,'_chromVAR_preproc_',chromVAR_preproc,'.RData',sep='')
  prmpt_var = paste('variability_chromVAR_norm_',chromVAR_norm,'_FIMO_',FIMO,'_same_',same,'_withPairs_',withPairs,'_binary_',binary,'_addBias_',bias_correct,'_chromVAR_preproc_',chromVAR_preproc,'.RData',sep='')
  
  source('compute_deviations_norm.R')
  save(dev,file=paste(chromVAR_org_pth,prmpt_dev,sep=''))
  if(bias_correct){
    variability = computeVariability(dev)
    save(variability,file=paste(chromVAR_org_pth,prmpt_var,sep=''))
  }
 
  print('Deviation Step Done!')
  
  if(!bias_correct){
    chromVAR_postproc=FALSE
  }
  
  prmpt = paste('cell_motif_chromVAR_postproc_',chromVAR_postproc,'_chromVAR_norm_',chromVAR_norm,'_FIMO_',FIMO,'_same_',same,'_withPairs_',withPairs,'_binary_',binary,'_addBias_',bias_correct,'_chromVAR_preproc_',chromVAR_preproc,'.RData',sep='')
  
  source('postproc.R')
  save(cell_motif,file=paste(chromVAR_org_pth,prmpt,sep=''))
  write.csv(cell_motif,file=paste(res_pth,'cell_motif.csv',sep=''),quote=FALSE,row.names=FALSE)
  mots = data.frame(motifs=colnames(cell_motif))
  write.csv(mots,file=paste(res_pth,'motifs.csv',sep=''),quote=FALSE,row.names=FALSE)  

  print('Postprocessing Step Done!')
  
  if (topic_modeling | chromVAR_downstream){
      if(chromVAR_downstream){
          source('chromVAR_downstream.R')
      }
      if(topic_modeling){
          source('cisTopic_preproc.R')
          source('create_cisTopicObject.R')
          source('run_cisTopic.R')
      }
      print('Downstream Step Done!')
  }
}
