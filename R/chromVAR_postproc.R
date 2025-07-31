source('remove_correlated_helper.R')
threshold = 0.90
vars <- variability[,'variability']
ix <- which(vars >= threshold)

ix2 <- ix[remove_correlated_helper(deviations(dev)[ix, , drop = FALSE],vars[ix])]
cell_motif <- t(deviations(dev)[ix2, , drop = FALSE])

