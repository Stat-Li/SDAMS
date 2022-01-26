
SDA <- function(sumExp, VOI = NULL, cleanData = TRUE, correctPValues = TRUE, ...){

    if (!is(sumExp, "SummarizedExperiment"))
        stop("Input must be an object of SummarizedExperiment class!")

    rawfeature <- assay(sumExp) 
    rawcoldata <- colData(sumExp)
    #newdata <- check_zero_data(rawfeature)
    if (cleanData) {
      newdata <- data_clean(rawfeature)
    } else {
      newdata <- list(feature = rawfeature,
                      feat.names = rownames(rawfeature))
    }
    newfeature <- newdata$feature
    
    newcoldata <- as(rawcoldata,"data.frame")
    feat.names <- newdata$feat.names

    rawresult <- lapply(seq_len(dim(newfeature)[1]), function (i) SDA.unit(
        featurevec=newfeature[i,], phenodata = newcoldata, VOI = VOI))
    results <- Reduce('comb', rawresult)
    if (correctPValues) {
      qv_1part <- apply(results$X1pvalue, 2, qvalue::qvalue,...)
      qv_2part <- qvalue::qvalue(results$X2pvalue,...)
    } else {
      qv_1part <- vector("list", 2)
      qv_1part[[1]] <- list(qvalues = rep(1, nrow(newfeature)))
      qv_1part[[2]] <- list(qvalues = rep(1, nrow(newfeature)))
      qv_2part <- list(qvalues = rep(1, nrow(newfeature)))
    }
    nparams = dim(results$pointest)[2]/2
    df.results <- list(gamma = as.matrix(results$pointest[,1:nparams]),
                        beta = as.matrix(results$pointest[,-(1:nparams)]),
                        pv_gamma = as.matrix(results$X1pvalue[,1]),
                        pv_beta = as.matrix(results$X1pvalue[,2]),
                        qv_gamma = as.matrix(qv_1part[[1]]$qvalues),
                        qv_beta = as.matrix(qv_1part[[2]]$qvalues),
                        pv_2part = results$X2pvalue,
                        qv_2part = qv_2part$qvalues,
                        feat.names = feat.names)
    return(df.results)
}

