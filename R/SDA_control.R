
SDA <- function(sumExp){

    if (!is(sumExp, "SummarizedExperiment"))
        stop("Input must be an object of SummarizedExperiment class!")

    rawfeature <- assay(sumExp); rawgroup <- colData(sumExp)$grouping
    newdata <- data_clean(rawfeature, rawgroup)
    newfeature <- newdata$feature
    newgroup <- newdata$group
    feat.names <- newdata$feat.names

    rawresult <- lapply(seq_len(dim(newfeature)[1]), function (i) SDA.unit(
        featurevec=newfeature[i,], grouping=newgroup))
    results <- Reduce('comb', rawresult)
    qv_1part <- apply(results$X1pvalue, 2, qvalue::qvalue)
    qv_2part <- qvalue::qvalue(results$X2pvalue)
    df.results <- list(gamma = results$pointest[,1],
                        beta = results$pointest[,2],
                        pv_gamma = results$X1pvalue[,1],
                        pv_beta = results$X1pvalue[,2],
                        qv_gamma = qv_1part[[1]]$qvalues,
                        qv_beta = qv_1part[[2]]$qvalues,
                        pv_2part = results$X2pvalue,
                        qv_2part = qv_2part$qvalues,
                        feat.names = feat.names)
    return(df.results)
}

