meme_motifs_pth = paste(workspace_dir_pth,'meme_motifs/',tolower(org),'.meme',sep='')
regions_fa_pth = paste(data_dir_pth,'regions.fa',sep='')
fimo_out_pth = paste(chromVAR_org_pth,'fimo_out/',sep='')

if (!('fimo_out'%in%list.files(chromVAR_org_pth))){
  system(paste('fimo -o',fimo_out_pth,'--max-stored-scores 2000000000',meme_motifs_pth,regions_fa_pth))
}
motif_lines = readLines(meme_motifs_pth)
motif_list = strsplit(motif_lines[grep('MOTIF ',motif_lines)],' ')
motifs = rep(NA,length(motif_list))
for (i in 1:length(motif_list)){
  motifs[i] = motif_list[[i]][3]
}
bed_file = read.table(bed_file_pth)
regions = paste(bed_file$V1,':',bed_file$V2,'-',bed_file$V3,sep='')
matches_pth = paste(fimo_out_pth,'fimo.tsv',sep='')
matches = read.table(matches_pth,sep='\t')
print('FIMO file read!')
colnames(matches) = matches[1,]
matches = matches[2:nrow(matches),]
fimo_bed_file = matches[,3:5]
fimo_bed_file_pth = paste(fimo_out_pth,'fimo_bed.bed',sep='')
write.table(fimo_bed_file,file=fimo_bed_file_pth,sep='\t',
            quote=FALSE,row.names=FALSE,col.names=FALSE)
fimo_intersect_bed_file_pth = paste(fimo_out_pth,'fimo_intersect.bed',sep='')

system(paste('bedtools intersect -wa -wb -a',fimo_bed_file_pth,'-b',bed_file_pth, '>',fimo_intersect_bed_file_pth))
system(paste('rm',fimo_bed_file_pth))
fimo_intersect_bed_file_pth = paste(fimo_out_pth,'fimo_intersect.bed',sep='')
intersect = read.table(fimo_intersect_bed_file_pth,sep='\t')
print('FIMO intersect file read!')
system(paste('rm',fimo_intersect_bed_file_pth))
seqnames = paste(intersect$V4,':',intersect$V5,'-',intersect$V6,sep='')
matches$seqnames = seqnames

cores = n_cores
cl <- makeCluster(cores) #not to overload your computer
registerDoParallel(cl)

partitions = round(seq(1, nrow(matches), by = nrow(matches)/n_seeds_cisTopic))[-1]
partitions = c(0,partitions,nrow(matches))


m_p <- foreach(i=2:(n_seeds_cisTopic+1),.combine='+') %dopar% {
  mat = matrix(0,nrow=length(regions),ncol=length(motifs))
  rownames(mat) = regions
  colnames(mat) = motifs
  for (j in (partitions[i-1]+1):partitions[i]){
    mot = matches[j,'motif_alt_id']
    seq = matches[j,'seqnames']
    if(mat[seq,mot]==0){
      mat[seq,mot] = 1
    }
    else{
      mat[seq,mot] = (mat[seq,mot] + 1)
    }
  }
  mat
  #do other things if you want
}
stopCluster(cl)

m_p = m_p[rownames(fragment_counts),]
m_p = as(m_p,'sparseMatrix')
m_p = m_p[,which(colSums(m_p)!=0)]
m_p_counts = m_p
save(m_p_counts,file=paste(chromVAR_org_pth,'motif_presence_FIMO_counts_dgC.RData',sep=''))
rm(m_p_counts)
m_p[which(m_p>1)] <- 1
save(m_p,file=paste(chromVAR_org_pth,'motif_presence_FIMO_dgC.RData',sep=''))
m_p = as(m_p,'lMatrix')
save(m_p,file=paste(chromVAR_org_pth,'motif_presence_FIMO_lgC.RData',sep=''))

source('convert_motif_ix_chromVAR.R')
rm(m_p)
