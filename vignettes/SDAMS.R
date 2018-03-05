### R code from vignette source 'SDAMS.Rnw'

###################################################
### code chunk number 1: quick start (eval = FALSE)
###################################################
# library("SDAMS")
# data("exampleSumExpset")
# results=SDA(exampleSumExpset)


###################################################
### code chunk number 2: directory (eval = FALSE)
###################################################
# path1 <- "/path/to/your/feature.csv/"
# path2 <- "/path/to/your/group.csv/"


###################################################
### code chunk number 3: GetDirectory
###################################################
directory1 <- system.file("extdata", package="SDAMS", mustWork=TRUE)
path1<-paste(directory1,"ProstateFeature.csv",sep="/")
directory2 <- system.file("extdata", package="SDAMS", mustWork=TRUE)
path2<-paste(directory2,"ProstateGroup.csv",sep="/")

###################################################
### code chunk number 4: CsvInput
###################################################
library("SDAMS")
exampleSEset1 = createSEsetFromCSV(path1,path2)
exampleSEset1

###################################################
### code chunk number 5: Accessors
###################################################
head(assay(exampleSEset1)[,1:10])
head(colData(exampleSEset1)$grouping)

###################################################
### code chunk number 6: MatrixInput
###################################################
set.seed(100)
featureInfo = matrix(runif(800,-2,5),ncol = 40)
featureInfo[featureInfo<0] = 0
rownames(featureInfo) = paste("feature",1:20,sep = '')
colnames(featureInfo) = paste('subject',1:40,sep = '')
groupInfo = data.frame(grouping=matrix(sample(0:1,40,replace = TRUE),ncol = 1))
rownames(groupInfo)=colnames(featureInfo)
exampleSEset2 = createSEsetFromEnvir(feature = featureInfo,group = groupInfo)
exampleSEset2
head(assay(exampleSEset2)[,1:10])
head(colData(exampleSEset2)$grouping)


###################################################
### code chunk number 7: results
###################################################
results = SDA(exampleSEset1)
head(results$gamma)
head(results$beta)
head(results$qv_gamma)
head(results$qv_beta)
head(results$qv_2part)
head(results$feat.names)


###################################################
### code chunk number 8: sessionInfo (eval = True)
###################################################
toLatex(sessionInfo())


