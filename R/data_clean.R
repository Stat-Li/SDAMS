
data_clean = function(rawfeature,rawgroup){
  ind = which(rowSums(rawfeature>0)>=10)
  feature = rawfeature[ind,]; group = as.matrix(rawgroup)
  feat.names = rownames(feature)
  return(list(feature=feature,group=group,feat.names=feat.names))
}

