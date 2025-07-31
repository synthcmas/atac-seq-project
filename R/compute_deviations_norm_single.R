compute_deviations_norm_single <- function(peak_set,
                                           counts_mat,
                                           background_peaks = NULL,
                                           expectation = NULL,
                                           fragments_per_sample = NULL,
                                           chromVAR_norm = TRUE,
                                           bias_correct = TRUE,
                                           intermediate_results = FALSE,
                                           threshold = 1) {
  
  if (length(peak_set) == 0) {
    return(list(z = rep(NA, ncol(counts_mat)), dev = rep(NA, ncol(counts_mat)),
                matches = 0, overlap = NA))
  }
  
  ### counts_mat should already be normed!
  tf_count <- length(peak_set)
  
  if (tf_count == 1) {
    observed <- as.vector(counts_mat[peak_set, ])
    
    if (chromVAR_norm){
      expected <- expectation[peak_set] * fragments_per_sample
      observed_deviation <- (observed - expected)/expected
    }
    else{
      observed_norm <- (observed/fragments_per_sample)*10000
    }
    
    if(bias_correct){
      sampled <- counts_mat[background_peaks[peak_set, ], ]
      sampled_expected <- outer(expectation[background_peaks[peak_set, ]],
                                fragments_per_sample)
      sampled_deviation <- (sampled - sampled_expected)/sampled_expected
      bg_overlap <- sum(background_peaks[peak_set,] == peak_set[1]) / 
        ncol(background_peaks)
    }
  } else {
    tf_vec <- sparseMatrix(j = peak_set,
                           i = rep(1, tf_count),
                           x = 1,
                           dims = c(1,
                                    nrow(counts_mat)))
    
    observed <- as.vector(tf_vec %*% counts_mat)
    
    if (!chromVAR_norm){
      observed_norm <- (observed/fragments_per_sample)*10000
    }else{
      expected <- as.vector(tf_vec %*% expectation %*% fragments_per_sample)
      observed_deviation <- (observed - expected)/expected
    }
    
    if(bias_correct){
      niterations <- ncol(background_peaks)
      sample_mat <-
        sparseMatrix(j = as.vector(background_peaks[peak_set, 
                                                    seq_len(niterations)]),
                     i = rep(seq_len(niterations), each = tf_count),
                     x = 1,
                     dims = c(niterations,
                              nrow(counts_mat)))
      
      sampled <- as.matrix(sample_mat %*% counts_mat)
      sampled_expected <- as.matrix(sample_mat %*% expectation %*% 
                                      fragments_per_sample)
      sampled_deviation <- (sampled - sampled_expected)/sampled_expected
      
      bg_overlap <- mean(sample_mat %*% t(tf_vec)) / tf_count      
    }
  }
  
  if(bias_correct){
    mean_sampled_deviation <- colMeans(sampled_deviation)
    sd_sampled_deviation <- apply(sampled_deviation, 2, sd)
    
    normdev <- (observed_deviation - mean_sampled_deviation)
    z <- normdev/sd_sampled_deviation
  }
  #fail_filter <- which(expected < threshold)
  

  
  #if (length(fail_filter) > 0) {
    #z[fail_filter] <- NA
    #normdev[fail_filter] <- NA
  #}
  
  if (intermediate_results) {
    if(bias_correct){
      out <- list(z = z, dev = normdev,
                  observed = observed,
                  sampled = sampled,
                  expected = expected,
                  sampled_expected = sampled_expected,
                  observed_deviation = observed_deviation,
                  sampled_deviation = sampled_deviation,
                  matches = tf_count / nrow(counts_mat),
                  overlap = bg_overlap)      
    }else{
      if(chromVAR_norm){
        out <- list(z = observed_deviation,
                    dev = observed_deviation,
                    observed = observed,
                    expected = expected,
                    observed_deviation = observed_deviation,
                    matches = tf_count / nrow(counts_mat),
                    overlap = NA)           
      }else{
        out <- list(z = observed_norm,
                    dev = observed_norm, 
                    observed = observed,
                    matches = tf_count / nrow(counts_mat),
                    overlap = NA)   
      }
    }
  } else {
    if(bias_correct){
      out <- list(z = z, dev = normdev,
                  matches = tf_count / nrow(counts_mat),
                  overlap = bg_overlap)      
    }else{
      if(chromVAR_norm){
        out <- list(z = observed_deviation,
                    dev = observed_deviation,
                    matches = tf_count / nrow(counts_mat),
                    overlap = NA)           
      }else{
        out <- list(z = observed_norm,
                    dev = observed_norm,
                    matches = tf_count / nrow(counts_mat),
                    overlap = NA)   
      }
    }
  }
  return(out)
}
