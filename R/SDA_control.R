
SDA <- function(sumExp, VOI = NULL, ...){

    if (!is(sumExp, "SummarizedExperiment"))
        stop("Input must be an object of SummarizedExperiment class!")

    rawfeature <- assay(sumExp); rawcoldata <- colData(sumExp)
    newdata <- data_clean(rawfeature)
    newfeature <- newdata$feature
    newcoldata <- as(rawcoldata,"data.frame")
    feat.names <- newdata$feat.names

    rawresult <- lapply(seq_len(dim(newfeature)[1]), function (i) SDA.unit(
        featurevec=newfeature[i,], phenodata = newcoldata, VOI = VOI))
    results <- Reduce('comb', rawresult)
    qv_1part <- apply(results$X1pvalue, 2, qvalue::qvalue,...)
    qv_2part <- qvalue::qvalue(results$X2pvalue,...)
    nparams <- dim(results$pointest)[2]/2
    df.results <- list(gamma = results$pointest[,1:nparams],
                        beta = results$pointest[,-(1:nparams)],
                        pv_gamma = results$X1pvalue[,1],
                        pv_beta = results$X1pvalue[,2],
                        qv_gamma = qv_1part[[1]]$qvalues,
                        qv_beta = qv_1part[[2]]$qvalues,
                        pv_2part = results$X2pvalue,
                        qv_2part = qv_2part$qvalues,
                        feat.names = feat.names)
    return(df.results)
    return(results)
}

