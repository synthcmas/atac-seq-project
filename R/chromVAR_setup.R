bed_file = read.table(bed_file_pth,sep='\t')
peaks = GRanges(seqnames = bed_file$V1,ranges = IRanges(start = bed_file$V2,
end = bed_file$V3))
peaks@ranges@start = peaks@ranges@start+as(1,'integer')
peaks@ranges@width = peaks@ranges@width-as(1,'integer')
# peaks = getPeaks(bed_file_pth,sort_peaks=TRUE)
rownames(mat) = NULL
fragment_counts <- SummarizedExperiment(assays = list(counts = mat),
                                        rowRanges = peaks)
peaks@ranges@start = peaks@ranges@start-as(1,'integer')
peaks@ranges@width = peaks@ranges@width+as(1,'integer')
rownames(fragment_counts) = as.character(peaks)
rm(peaks)
rm(mat)
colData(fragment_counts) <- DataFrame(cell_types = labels)
fragment_counts <- filterPeaks(fragment_counts, non_overlapping = TRUE)
