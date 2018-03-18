
######## read data (two ways)   ########
#-------- get data from seperate matrix ---------
createSEFromMatrix <- function(feature, group) {

    subjectName = rownames(group)
    if (is(group, "numeric")) {
        group <- as.data.frame(group)
        colData <- as(group,'DataFrame')
        rownames(colData) <- subjectName
        colnames(colData)[1] <- "grouping"
    }
    if (is(group, "data.frame")) {
        mgroup <- model.matrix(~group[,1])
        group <- as.data.frame(mgroup[,2])
        colData <- DataFrame(grouping=group)
        rownames(colData) <- subjectName
        colnames(colData)[1] <- "grouping"
    }
    if (is(group, "matrix")) {
        group <- as.data.frame(group)
        colData <- as(group,'DataFrame')
        rownames(colData) <- subjectName
        colnames(colData)[1] <- "grouping"
    }

    feature <- as.matrix(feature)
    if (ncol(feature)!=nrow(group)) {
        stop("Feature data and group data do not match!")
    }


    result <- SummarizedExperiment(assays = SimpleList(counts=feature),
                                   colData = colData)
    return(result)
}

# ------- import data from csv files ----------
createSEFromCSV <- function(featurePath, groupPath, rownames1 = 1,
                            rownames2 = 1, header1 = TRUE, header2 = TRUE){

    feature <- read.csv(featurePath, row.names = rownames1, header = header1,
                        check.names = FALSE)
    group <- read.csv(groupPath, row.names = rownames2, header = header2,
                        check.names = FALSE)

    result <- createSEFromMatrix(feature = feature, group = group)

    return(result)
}
