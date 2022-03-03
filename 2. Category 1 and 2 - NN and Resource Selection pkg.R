#Loading required packages
library(nloptr)
library(msm)
library(ResourceSelection)
library(nnet)

#####################################################################
##LK Method
###The model fitting for LK method adopted from the R code 
#provided by Phillips and Elith (2013)
LK <- function(presence, background) {
  lik <- function(beta) {
    datapsi <- plogis(presence %*% beta)
    gridpsi <- plogis(background %*% beta)
    -mean(log(datapsi)) + log(sum(gridpsi))
  }
  presence <- cbind(1, as.matrix(presence))
  background <- cbind(1, as.matrix(background))
  nlm(lik, rep(0, NCOL(presence)))$estimate
}
########################################################################
#CLK Method
NEW<-function(prevalence,presence,background,fun,quadratic,cubic)
{
  if(cubic){
    npar <- 4
  } 
  else if (quadratic) {
    npar <- 3
  }
  else {
    npar <- 2
  }
  presence <- cbind(1, as.matrix(presence))
  background <- cbind(1, as.matrix(background))
  x0<-c(rep(1,npar))
  fn<-function(x) -mean(log(fun(presence,x)))
  heq <- function(x) mean(fun(background,x))-prevalence # heq == 0
  auglag(c(rep(0.5,npar)), fn, heq = heq,localsolver="lbfgs")$par 
  # with COBYLA  choose c(1,1) as the starting value
} 
#####################################################################

fun0<-function (x,beta) {plogis(x%*%beta)}  ##logit function
fun1 <- function(x,beta){exp(x%*%beta)} #exponential

#####################################################################

#Function for plotting probability lines
makeline <- function(beta, dat,fun, ...)
{
  preds <- apply(as.matrix(dat), 1, function(x) fun(c(1,x),beta))
  lines(as.matrix(dat)[,1], preds, ...)
}

#####################################################################
#Calculate the average of the betas
mean_preds <- function (beta, dat,...){
  preds <- apply(as.matrix(dat), 1, function(x) plogis(beta %*% c(1,x)))
  preds_data<-preds[sapply(preds, is.numeric)]  
  sapply(preds_data, mean, na.rm = T)
}

#######################################################################
# It takes two required arguments:
#   "pp" is the probability of presence as a function of one predictor
#   "alg" is a function that calculates the parameters of a logistic model
# The "alg" function must take two parameters: a matrix of covariates at 
# presence sites and a matrix of covariates at background sites.
# Optional arguments are:
#   "quadratic" -- true if quadratic terms are to be used
#   "nlines"  -- the number of replicates in the experiment
#   "xl" -- true if the x axis labels are to be included in the plot
#  "x1" and "x2" limits for the range of interest to be used/ the 
#   range in which the function is. 
#Has been adapted from Phillips and Elith (2013)

makeplot <- function(pp, fun,x1, x2, quadratic=F, cubic=F,nlines=1000) {
  np <- 5000 # presence samples
  nb <- 50000 #Background samples
  
  # generate test data
  testx <- sort(seq(x1, x2, length.out=1000) )
  
  truth <- pp(testx)#Original Function
  
  # make some presence-absence data
  y1 <- rbinom(1000, 1, truth)
  
  if (cubic){
    test <- data.frame(x=testx,xx=testx*testx,xxx=testx*testx*testx)
  }
  else if (quadratic){
    test <- data.frame(x=testx, xx=testx*testx)
  }
  else{
    test <- data.frame(x=testx)
  }

  #Defining required data frames
  s_test <- y1
  ua_test <- test
  beta_list <- data.frame()
  CLK_list <- data.frame()
  
  for (line in 1:nlines) {
    np <- 5000       # number of presence points
    nb <- 50000      # number of background points
    nt <- 100000
    
    # generate presence data
    x <-  seq(x1, x2, length.out = nb)
    p <- pp(x)
    
    # make some presence-absence data
    y <- rbinom(nt, 1, p)
    if (sum(y) < np) stop("Prevalence is low, need more than 50000 samples") 
    
    # keep np presence points
    px <- x[sample.int(nb,np,replace=TRUE,prob=pp(x))]
    if (cubic){
      presence <- data.frame(x=px,xx=px*px,xxx=px*px*px)
    }
    else if (quadratic){
      presence <- data.frame(x=px, xx=px*px)
    }
    else{
      presence <- data.frame(x=px)
    }
    # generate random background data
    bx <- seq(x1, x2, length.out=nb) 
    
    if (cubic){
      background <- data.frame(x=bx,xx=bx*bx,xxx=bx*bx*bx)
    }
    else if (quadratic){
      background <- data.frame(x=bx, xx=bx*bx)
    }
    else{
      background <- data.frame(x=bx)
    }
    #################
    #LK Method
    #Fitting using the ResourceSelection package
    
    ynew <- c(rep(1,NROW(presence)), rep(0,NROW(background)))
    xnew <- rbind(presence,background)
    data_LK <- data.frame(ynew,xnew)
    
    if (cubic){
      model_LK <- rspf(ynew~x + xx + xxx,data_LK,m=0,B=0)
    }
    else if (quadratic){
      model_LK <- rspf(ynew~x + xx,data_LK,m=0,B=0)
    }
    else{
      model_LK <- rspf(ynew~x,data_LK,m=0,B=0)
    }
    
    coef <- as.data.frame(model_LK$coefficients)
    
    if (cubic){
      y <- plogis(coef[1,] + coef[2,]*testx + coef[3,]*(testx*testx) + coef[4,]*(testx*testx*testx))
      beta_list <- rbind(beta_list,t(coef))
    }
    else if (quadratic){
      y <- plogis(coef[1,] + coef[2,]*testx + coef[3,]*(testx*testx))
      beta_list <- rbind(beta_list,t(coef))
    }
    else{
      y <- plogis(coef[1,] + coef[2,]*testx)
      beta_list <- rbind(beta_list,t(coef))
    }
    
    ######Using the LK Method function above 
    #(Optional-can change to this if needed)
    #if(cubic){
    # npar <- 4
    #} 
    #else if (quadratic) {
    #npar <- 3
    #}
    #else {
    # npar <- 2
    #}
    #LK_beta <-matrix(0,nlines,npar)
    
    #LK_beta[line,] <- LK(presence,background)
    #makeline(LK_beta[line,], test,fun0, col="#d9ef8b", lty = 2)
    
    #beta_list <- rbind(beta_list,t(LK_beta[line,]))
    ################
    #Estimating C (Using PBL Method)
    s <- c(rep(1,NROW(presence)), rep(0,NROW(background)))
    ua <- rbind(presence, background)
    
    #Fit model for the conditional probability of selection
    model_nnet<-nnet(s~., data=data.frame(cbind(s,ua)),size=5,range=0.1,decay=5e-2)
    nnet_pred<-predict(model_nnet)
    
    #estimate c using ranked predictions
    est_c<-mean(nnet_pred[order(-nnet_pred)][1:(nb*0.1)])
    
    #################
    #CLK Method
    if(cubic){
      npar <- 4
    } 
    else if (quadratic) {
      npar <- 3
    }
    else {
      npar <- 2
    }
    NEWbeta0<-matrix(0,nlines,npar)
    
    p0_mean <- np*(1-est_c)/est_c
    prevalence <- (p0_mean)/nb
    
    #Fit the CLK method
    NEWbeta0[line,]<-NEW(prevalence, presence, background,fun,quadratic, cubic)
    
    CLK_list <- rbind(CLK_list,t(NEWbeta0[line,]))
    
  }
  write.csv(beta_list, "LK_Beta_saved.csv")
  write.csv(CLK_list, "CLK_Beta_saved.csv")
}

#Source the file that has slected functions
source("2. Selected functions of the categories.R")

#####################################################################
#Logistic functions
#Category 1
pdf("fig - LC and RSPF satisfied Category 1.pdf", width=8, height=12)
op <- par(mfrow=c(1,3), las=1, cex=1.2, mar=c(3,4,4,2)+0.1)
makeplot(logistic_3 , x1=-1,x2=1,quadratic=F, cubic=F, fun=fun0)
makeplot(Qlogistic_3 , x1=-1,x2=1, quadratic=T, cubic=F, fun=fun0)
makeplot(Clogistic_4 , x1=-1,x2=1, quadratic=F, cubic=T, fun=fun0)
par(op)
dev.off()

#####################################################################
#Category 2
pdf("fig - Only LC Satisifed Category 2.pdf", width=8, height=12)
op <- par(mfrow=c(1,3), las=1, cex=1.2, mar=c(3,4,4,2)+0.1)
makeplot(semilogistic, fun = fun1, x1=0,x2=1, quadratic = F, cubic = F)
makeplot(gaussian, fun = fun1, x1=0,x2=1, quadratic = T, cubic = F)
makeplot(exponent, fun = fun1, x1=0,x2=1, quadratic = F, cubic = T)
par(op)
dev.off()

#####################################################################
#####################################################################

#Cloglog Functions 
#Category 1 (Appendix)
makeplot(cloglog_1, x1=-1,x2=1, quadratic=F, cubic=F, fun=fun0)
makeplot(Qcloglog_1 , x1=-1,x2=1, quadratic=T, cubic=F, fun=fun0)
makeplot(Ccloglog_2 , quadratic=F, cubic=T, fun=fun0)

