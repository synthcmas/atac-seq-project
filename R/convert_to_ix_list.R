convert_to_ix_list <- function(ix) {
  stopifnot(inherits(ix, "Matrix") || inherits(ix, "matrix"))
  if (inherits(ix, "sparseMatrix")) {
    tmp <- summary(ix)
    tmp$j <- factor(tmp$j, levels = seq_len(ncol(ix)), ordered = TRUE)
    out <- split(tmp$i, tmp$j)
  } else {
    out <- lapply(seq_len(ncol(ix)), function(x) which(ix[, x]))
  }
  names(out) <- colnames(ix)
  return(out)
}