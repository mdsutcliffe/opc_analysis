runPagoda <- function(rsem, # rsem counts, including annotation
                      iControls, # indices of split-pool controls
                      iSamples # indices of ten-cell samples
                      ) {

  samples <- rsem[,iSamples]
  controls <- rsem[,iControls]
  
  samples <- data.frame(samples)
  controls <- data.frame(controls)
  
  knn_samples <- knn.error.models(samples,
                                  k = ncol(samples)/4,
                                  n.cores = 1,
                                  min.count.threshold = 5,
                                  min.nonfailed = 5,
                                  save.model.plots = F)
  knn_controls <- knn.error.models(controls,
                                   k = ncol(controls)/4,
                                   n.cores = 1,
                                   min.count.threshold = 5,
                                   min.nonfailed = 5,
                                   save.model.plots = F)
  
  varinfo_samples <- pagoda.varnorm(knn_samples, counts = samples, trim = 3/ncol(samples), n.cores=1, plot = F)
  varinfo_controls <- pagoda.varnorm(knn_controls, counts = controls, trim = 3/ncol(controls),n.cores=1, plot = F)
  
  # Control for sequencing depth
  varinfo_samples <- pagoda.subtract.aspect(varinfo_samples, colSums(samples[, rownames(knn_samples)]>0))
  varinfo_controls <- pagoda.subtract.aspect(varinfo_controls, colSums(controls[, rownames(knn_controls)]>0))
  
  # Evaluate distributions of adj_var for both controls and samples
  samples_var <- data.matrix(varinfo_samples$arv)
  colnames(samples_var) <- "samples_var"
  controls_var <- data.matrix(varinfo_controls$arv)
  colnames(controls_var) <- "controls_var"
  merge_var <- merge(samples_var,controls_var,by='row.names')
  
  thresh <- quantile(merge_var$controls_var, 0.95)
  
  candidate_genes <- as.character(merge_var[merge_var$samples_var > thresh & merge_var$controls_var <= thresh,1])
  
  merge_var$candidate <- merge_var[,1] %in% candidate_genes
  
  return(list("genes" = candidate_genes, 
              "adjusted_variance" = merge_var,
              "varinfo_controls" = varinfo_controls,
              "varinfo_samples" = varinfo_samples,
              "thresh" = thresh))
}
