remove_correlated_helper <- function(mat, val, cutoff = 0.90) {
  stopifnot(nrow(mat) == length(val))
  cormat <- cor(t(mat), use = "pairwise.complete.obs")
  diag(cormat) <- NA
  keep <- seq_len(nrow(mat))
  for (i in order(val, decreasing = TRUE)) {
    if (i %in% keep) {
      toremove <- which(cormat[keep, i] >= cutoff)
      if (length(toremove) > 0) 
        keep <- keep[-toremove]
    }
  }
  return(keep)
}
