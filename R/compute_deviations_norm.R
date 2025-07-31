source('convert_to_ix_list.R')
source('compute_deviations_norm_single.R')
register(MulticoreParam(n_cores))

if (bias_correct){
  set.seed(123)
  background_peaks = getBackgroundPeaks(fragment_counts)
}else{
  background_peaks = NULL
}

coldata = colData(fragment_counts)
rowData = colData(motif_ix)

peak_indices = convert_to_ix_list(annotationMatches(motif_ix))
rm(motif_ix)
counts_mat = counts(fragment_counts)
expectation = computeExpectations(fragment_counts)
rm(fragment_counts)
fragments_per_sample = colSums(counts_mat)
sample_names <- colnames(counts_mat)

results <- bplapply(peak_indices,
                    compute_deviations_norm_single,
                    counts_mat = counts_mat,
                    background_peaks = background_peaks,
                    chromVAR_norm = chromVAR_norm,
                    bias_correct = bias_correct,
                    expectation = expectation,
                    fragments_per_sample = fragments_per_sample)

rm(fragments_per_sample)
rm(expectation)
rm(peak_indices)
rm(background_peaks)

z <- t(vapply(results, function(x) x[["z"]], rep(0, ncol(counts_mat))))
dev <- t(vapply(results, function(x) x[["dev"]], rep(0, ncol(counts_mat))))
rm(counts_mat)
colnames(z) <- colnames(dev) <- sample_names
rm(sample_names)
rowData$fractionMatches <- vapply(results, function(x) x[["matches"]], 0)
rowData$fractionBackgroundOverlap <- vapply(results,
                                            function(x) x[["overlap"]],
                                            0)
out <- SummarizedExperiment(assays = list(deviations = dev, z = z),
                            colData = coldata,
                            rowData = rowData)
rm(z)
dev = new("chromVARDeviations", out)
rm(out)
