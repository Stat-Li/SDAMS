## generate feature and phenotype group data
featureInfo = matrix(runif(800,-2,5),ncol = 40)
featureInfo[featureInfo<0] = 0
feat.names = paste("feature",1:20,sep = '')
rownames(featureInfo) = feat.names
sub.names = paste('subject',1:40,sep = '')
colnames(featureInfo) = sub.names
groupInfo = data.frame(grouping=matrix(sample(0:1,40,replace = TRUE),ncol = 1))
rownames(groupInfo)=colnames(featureInfo)

## input data to create a "MSset" object

exampleMSset2 = createMSsetFromEnvir(feature = featureInfo,group = groupInfo)

expect_s4_class(exampleMSset2,'MSset')
expect_true(all(featuredata(exampleMSset2) == featureInfo))
expect_true(all(phenotypedata(exampleMSset2) == groupInfo))
expect_true(all(rownames(featuredata(exampleMSset2)) == feat.names))
expect_true(all(colnames(featuredata(exampleMSset2)) == sub.names))







