
######## read data (two ways)   ########
#-------- get data from seperate matrix ---------
createSEFromMatrix <- function(feature, colData) {

    subjectName = rownames(colData)
    cName = colnames(colData)

    if (is(colData, "data.frame")) {
        colData <- as(colData, "DataFrame")
    }
    if (is(colData, "matrix")) {
        colData <- as(colData,'DataFrame')
    }

    feature <- as.matrix(feature)
    if (ncol(feature)!=nrow(colData)) {
        stop("Feature data and column data do not match!")
    }


    result <- SummarizedExperiment(assays = SimpleList(counts=feature),
                                   colData = colData)
    return(result)
}

# ------- import data from csv files ----------
createSEFromCSV <- function(featurePath, colDataPath, rownames1 = 1,
                            rownames2 = 1, header1 = TRUE, header2 = TRUE){

    feature <- read.csv(featurePath, row.names = rownames1, header = header1,
                        check.names = FALSE)
    colData <- read.csv(colDataPath, row.names = rownames2, header = header2,
                        check.names = FALSE)
    
    result <- createSEFromMatrix(feature = feature, colData = colData)

    return(result)

}
