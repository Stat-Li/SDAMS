

comb <- function(list1, list2){
    res <- list()

    res$pointest <- rbind(list1$pointest, list2$pointest)
    res$X1pvalue <- rbind(list1$X1pvalue, list2$X1pvalue)
    res$X2pvalue <- rbind(list1$X2pvalue, list2$X2pvalue)

    return(res)
}
