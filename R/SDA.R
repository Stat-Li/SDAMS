
aft_model1<- function(data_feature, data_group,binit=0,bw=1){
  x <- data_group; y <- data_feature
  n <- length(y)

  if (bw==1){an <- 1.144*sd(log(y)-as.vector(x%*%as.matrix(binit)))*n^(-1/5)}
  if (bw==2){an <- sd(log(y)-as.vector(x%*%as.matrix(binit)))*n^(-1/5)}
  if (bw==3){an <- sd(log(y)-as.vector(x%*%as.matrix(binit)))*n^(-1/7)}
  if (bw==4){an <- sd(log(y)-as.vector(x%*%as.matrix(binit)))*n^(-1/9)}
  if (bw==5){an <- 4^(1/3)*min(sd(log(y)-as.vector(x%*%as.matrix(binit))),
                              IQR(log(y)-as.vector(x%*%as.matrix(binit)))/1.34)*n^(-1/5)}
  if (bw==6){an <- (8*sqrt(2)/3)^(1/5)*min(sd(log(y)-as.vector(x%*%as.matrix(binit))),
                                          IQR(log(y)-as.vector(x%*%as.matrix(binit)))/1.34)*n^(-1/5)}
  if (bw==7){an <- 4^(1/3)*min(sd(log(y)-as.vector(x%*%as.matrix(binit))),
                              IQR(log(y)-as.vector(x%*%as.matrix(binit)))/1.34)*n^(-1/3)}

  kern <- dnorm ##kernel function: standard normal function
  kern.1st <- function(x){-x*dnorm(x)}
  kern.2nd <- function(x){(x^2-1)*dnorm(x)}
  kern.cdf <- pnorm

  ##-------------paired difference between each error----------------------
  e_diff <- function(beta.iter,x,y){
    e.diff<- outer(-log(y)+as.vector(t(x%*%as.matrix(beta.iter))),
                   log(y)-as.vector(t(x%*%as.matrix(beta.iter))), '+')
    return(e.diff)
  }

  loglikf_1 <- function(beta.iter, x, y){
    e.diff<- e_diff(beta.iter,x,y)
    loglik_value<- -1*sum(log(y))+sum(log(rowSums(kern(e.diff/an))/n/an))
    return(loglik_value)
  }

  fbeta_1 <- function(beta.iter, x, y){
    e.diff<- e_diff(beta.iter,x,y)
    fbeta<- rep(NA, length(beta.iter))
    for (p in 1:length(beta.iter)) {
      x_vec<- as.matrix(x)[,p]
      x.diff<- outer(x_vec,-x_vec,'+')
      fbeta[p]<- 1/an*sum(rowSums(kern.1st(e.diff/an)*x.diff)/rowSums(kern(e.diff/an)))
    }
    return(fbeta)
  }

  fbeta_dev_1 <- function(beta.iter,x,y){
    e.diff <- e_diff(beta.iter,x,y)
    hessian.m <- matrix(NA, nrow = length(beta.iter), ncol = length(beta.iter))
    for (i in 1:length(beta.iter)) {
      x_i<- as.matrix(x)[,i]
      x_i_diff<- outer(x_i,-x_i,'+')
      for (j in 1:length(beta.iter)) {
        x_j<- as.matrix(x)[,j]
        x_j_diff<- outer(x_j,-x_j,'+')
        hessian.m[i,j]<- 1/an/an*sum((rowSums(kern.2nd(e.diff/an)*x_i_diff*x_j_diff)*rowSums(kern(e.diff/an))-
                                        rowSums(kern.1st(e.diff/an)*x_i_diff)*rowSums(kern.1st(e.diff/an)*x_j_diff))/(rowSums(kern(e.diff/an)))^2)
      }
    }
    return(hessian.m)
  }

  objfun_1 <- function(beta.iter){
    stopifnot(is.numeric(beta.iter))
    f <- loglikf_1(beta.iter,x,y)
    g <- fbeta_1(beta.iter,x,y)
    B <- as.matrix(fbeta_dev_1(beta.iter,x,y))
    list(value =f, gradient = g, hessian = B)
  }

  trust.results <- try(trust(objfun_1, binit, 1, 5, minimize = FALSE),silent = TRUE)
  if (class(trust.results) == 'try-error') {
    return(list(pointest=NA,seest=NA,null.deviance=NA,residual.deviance=NA))
  }
  else {
    beta.est <- trust.results$argument
    se.est <- try(sqrt(-1*solve(trust.results$hessian))[1,1],silent = TRUE)
    if (class(se.est) == 'try-error'){se.est=NA;null.deviance=0;residual.deviance=0}
    else {
      se.est=se.est
      null.deviance = -2*loglikf_1(0, x, y)
      residual.deviance = -2*loglikf_1(beta.est, x, y)
    }
    return(list(pointest=beta.est,seest=se.est,null.deviance=null.deviance,
                residual.deviance=residual.deviance))
  }
}


SDA.unit = function(featurevec, grouping,bw=1){
  data0 = data.frame(featurevec, grouping)
  data_binary = data0; data_binary$featurevec[data_binary$featurevec>0] = 1
  data_AFT = data0[data0$featurevec>0,]

  non0_cnt <- c(sum(data_binary$featurevec[data_binary$grouping==0]==1),
                     sum(data_binary$featurevec[data_binary$grouping==1]==1))
  zero_cnt <- c(sum(data_binary$featurevec[data_binary$grouping==0]==0),
                  sum(data_binary$featurevec[data_binary$grouping==1]==0))
  group_size <- c(sum(data0$grouping==0),sum(data0$grouping==1))

  if(any(non0_cnt==group_size)){
    coef_logit=NA; se_logit=NA
    aft_summary <- aft_model1(data_AFT$featurevec, data_AFT$grouping,bw=bw)
    coef_aft <- aft_summary$pointest; se_aft <- aft_summary$seest
    diff.dev.logit = 0
    diff.dev.aft = aft_summary$null.deviance-aft_summary$residual.deviance
  }
  else if(any(non0_cnt<2) || all(data_AFT$featurevec[1]==data_AFT$featurevec)){
    coef_aft=NA; se_aft=NA
    logit_reg <- glm(data_binary$featurevec ~ data_binary$grouping, family = 'binomial', control = list(maxit = 50))
    logit_summary <- summary(logit_reg)
    coef_logit <- logit_summary$coefficients[2,1]; se_logit <- logit_summary$coefficients[2,2]
    diff.dev.logit = logit_summary$null.deviance-logit_summary$deviance
    diff.dev.aft = 0
  }
  else{
    logit_reg <- glm(data_binary$featurevec ~ data_binary$grouping, family = 'binomial', control = list(maxit = 50))
    logit_summary <- summary(logit_reg)
    coef_logit <- logit_summary$coefficients[2,1]; se_logit <- logit_summary$coefficients[2,2]
    aft_summary <- aft_model1(data_AFT$featurevec, data_AFT$grouping,bw=bw)
    coef_aft <- aft_summary$pointest; se_aft <- aft_summary$seest
    diff.dev.logit = logit_summary$null.deviance-logit_summary$deviance
    diff.dev.aft = aft_summary$null.deviance-aft_summary$residual.deviance
  }

  if (diff.dev.logit == 0) {p_logit = NA}
  else {p_logit = 1-pchisq(diff.dev.logit,1)}
  if (diff.dev.aft == 0) {p_aft = NA}
  else {p_aft = 1-pchisq(diff.dev.aft,1)}

  diff.total = diff.dev.logit + diff.dev.aft

  if ((diff.dev.logit==0)|(diff.dev.aft==0)) {
    X2pv = 1-pchisq(diff.total,1)
  } else {
    X2pv = 1-pchisq(diff.total,2)
  }

  return(list(pointest=c(coef_logit, coef_aft),X1pvalue=c(p_logit,p_aft),X2pvalue=X2pv))
}

