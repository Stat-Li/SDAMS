
data_clean <- function(rawfeature){

    ind <- which(rowSums(rawfeature>0)>=10)
    feature <- rawfeature[ind, ]
    feat.names <- rownames(feature)

    return(list(feature = feature,
                feat.names = feat.names)
            )
}

