#rm(peaks)
#rm(mat)
load(lib_size_pth)
colData(fragment_counts)$depth = lib_size
counts_mat = counts(fragment_counts)
min_in_peaks = 0.5*median(colSums(counts_mat)/lib_size)
rm(counts_mat)
min_depth = max(500,0.1*median(lib_size))
rm(lib_size)
fragment_counts <- filterSamples(fragment_counts, min_depth = min_depth, 
                                 min_in_peaks = min_in_peaks, shiny = FALSE)
labels = labels[colnames(fragment_counts)]
fragment_counts = filterPeaks(fragment_counts,non_overlapping=TRUE)
