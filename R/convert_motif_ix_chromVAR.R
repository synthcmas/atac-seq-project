save(m_p,file=paste(chromVAR_org_pth,'m_p.RData',sep=''))
coldata = data.frame(name=colnames(m_p))
rownames(coldata) = colnames(m_p)
matches = list()
matches[['motifMatches']] = m_p


motif_ix = SummarizedExperiment(assays=matches,
                                rowRanges=rowRanges(fragment_counts),
                                colData=coldata)
rowData(motif_ix)=rowData(fragment_counts)


