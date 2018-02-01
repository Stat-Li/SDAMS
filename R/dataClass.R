setClass("MSset",
         contains = "ExpressionSet")

# set up for "feature"
setGeneric("featuredata", function(object) standardGeneric("featuredata"))
setGeneric("featuredata<-", function(object, value) standardGeneric(
  "featuredata<-"))
setMethod("featuredata", signature(object = "MSset"), function(object) {
  exprs(object)
})
setReplaceMethod("featuredata", signature(object = "MSset", value = "matrix"),
                 function(object, value) {
                   exprs(object) <- value
                   object
                 })

# set up for "group": phenotype data
setGeneric("phenotypedata", function(object) standardGeneric("phenotypedata"))
setGeneric("phenotypedata<-", function(object, value) standardGeneric(
  "phenotypedata<-"))
setMethod("phenotypedata", signature(object = "MSset"), function(object) {
  pData(object)
})
setReplaceMethod("phenotypedata", signature(object = "MSset", value = "data.frame"),
                 function(object, value) {
                   pData(object) <- value
                   object
                 })
######## read data (two ways)   ########
#-------- get data from current environment ---------
createMSsetFromEnvir <- function(feature, group) {

  if (is(group, "numeric")) {
    group <- as.data.frame(group)
    colnames(group) <- 'grouping'
  }
  if (is(group, "data.frame")) {
    mgroup <- model.matrix(~group[,1])
    group <- as.data.frame(mgroup[,2])
    colnames(group) <- 'grouping'
  }
  if (is(group, "matrix")) {
    group <- as.data.frame(group)
    colnames(group) <- 'grouping'
  }

  feature <- as.matrix(feature)
  rownames(group) = colnames(feature)

  if (is(group, "data.frame") || is(group, "AnnotatedDataFrame")) {
    group <- as(group, "AnnotatedDataFrame")
  }
  result <- new("MSset", exprs = feature, phenoData = group)
  return(result)
}

# ------- import data from csv files ----------
createMSsetfromCSV <- function(featurePath, groupPath, rownames1=1,rownames2=1,
                            header1=TRUE,header2=TRUE){

  feature <- read.csv(featurePath, row.names = rownames1, header = header1,
                      check.names = FALSE)
  group <- read.csv(groupPath, row.names = rownames2, header = header2,
                    check.names = FALSE)

  if (is(group, "numeric")) {
    group <- as.data.frame(group)
    colnames(group) <- 'grouping'
  }
  if (is(group, "data.frame")) {
    mgroup <- model.matrix(~group[,1])
    group <- as.data.frame(mgroup[,2])
    colnames(group) <- 'grouping'
  }
  if (is(group, "matrix")) {
    group <- as.data.frame(group)
    colnames(group) <- 'grouping'
  }

  feature <- as.matrix(feature)
  rownames(group) = colnames(feature)

  if (is(group, "data.frame") || is(group, "AnnotatedDataFrame")) {
    group <- as(group, "AnnotatedDataFrame")
  }
  result <- new("MSset", exprs = feature, phenoData = group)
  return(result)
}
