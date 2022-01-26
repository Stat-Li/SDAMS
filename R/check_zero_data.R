
check_zero_data <- function(rawfeature){
  n_sample <- ncol(rawfeature)
  half_sample <- ceiling(n_sample / 4)
  low_ind <- rowSums(rawfeature == 0) > half_sample
  if (sum(low_ind) > 0) {
    warning_msg = paste0("There are ", sum(low_ind), " features of ", length(low_ind), " where more than 1/2 of the samples are zero.\nProceed with caution!", collapse = "")
    warning(warning_msg)
  }
  feature <- rawfeature[!low_ind, ]
  feat.names <- rownames(feature)

  return(list(feature = feature,
              feat.names = feat.names)
          )
}

data_clean <- function(rawfeature) {
  ind <- which(rowSums(rawfeature > 0) >= 10)
  feature <- rawfeature[ind, ]
  feat.names <- rownames(feature)
  return(list(feature = feature,
              feat.names = feat.names)
  )
  
}
