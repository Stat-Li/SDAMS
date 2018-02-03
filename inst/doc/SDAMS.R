### R code from vignette source 'SDAMS.Rnw'

###################################################
### code chunk number 1: quick start (eval = FALSE)
###################################################
# library("SDAMS")
# data("exampleMSset")
# results=SDA(exampleMSset)


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
exampleMSset1 = createMSsetfromCSV(path1,path2)
exampleMSset1
head(featuredata(exampleMSset1)[,1:10])
head(phenotypedata(exampleMSset1))


###################################################
### code chunk number 5: MatrixInput
###################################################
set.seed(100)
featureInfo = matrix(runif(800,-2,5),ncol = 40)
featureInfo[featureInfo<0] = 0
rownames(featureInfo) = paste("feature",1:20,sep = '')
colnames(featureInfo) = paste('subject',1:40,sep = '')
groupInfo = data.frame(grouping=matrix(sample(0:1,40,replace = TRUE),ncol = 1))
rownames(groupInfo)=colnames(featureInfo)
exampleMSset2 = createMSsetFromEnvir(feature = featureInfo,group = groupInfo)
exampleMSset2
head(featuredata(exampleMSset2)[,1:10])
head(phenotypedata(exampleMSset2))



###################################################
### code chunk number 6: result1
###################################################
results = SDA(exampleMSset1)
head(results$gamma)
head(results$beta)
head(results$qv_gamma)
head(results$qv_beta)
head(results$qv_2part)
head(results$feat.names)


###################################################
### code chunk number 7: sessionInfo (eval = True)
###################################################
toLatex(sessionInfo())


