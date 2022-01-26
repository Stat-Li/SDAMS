## generate a matrix with feature and  group data
set.seed(100)
featureInfo <- matrix(runif(8000, -2, 5), ncol = 40)
featureInfo[featureInfo<0] <- 0
rownames(featureInfo) <- paste("gene", 1:200, sep = '')
colnames(featureInfo) <- paste('cell', 1:40, sep = '')
n_zero <- rowSums(featureInfo == 0)
more_half <- which(n_zero > 20)

for (i_more in more_half) {
  which_zero <- which(featureInfo[i_more, ] == 0)
  choose_change <- sample(which_zero, length(which_zero) - 20)
  featureInfo[i_more, choose_change] <- runif(length(choose_change), 1, 5)
}


groupInfo <- data.frame(grouping=matrix(sample(0:1, 40, replace = TRUE),
                                        ncol = 1))
rownames(groupInfo) <- colnames(featureInfo)

exampleSE2 <- createSEFromMatrix(feature = featureInfo, colData = groupInfo)
exampleSE2List <- list(feature = featureInfo, groups = groupInfo)

test_that("using a list fails", {
  # using a list fails
  expect_error(SDA(exampleSE2List))
})

test_that("using the SE object works", {
  # using the object works
  expect_snapshot(SDA(exampleSE2))
})

test_that("having a feature with > 20 zeros causes warning", {
  featureInfo[1, sample(ncol(featureInfo), 21)] <- 0
  exampleSE3 <- createSEFromMatrix(feature = featureInfo, colData = groupInfo)
  expect_snapshot_warning(SDA(exampleSE3))
})

