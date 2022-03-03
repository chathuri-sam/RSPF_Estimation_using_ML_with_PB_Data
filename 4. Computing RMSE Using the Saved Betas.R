#This file needs the respective otputs from the code in
#3. Category 1 and 2 - NN and Resource Selection pkg and
# 4. Category 3 NN and Resource Selection pkg with c 
# for each of the functions in the categories.

#Plotting the probability lines
makeline <- function(beta, dat,fun, ...)
{
  preds <- apply(as.matrix(dat), 1, function(x) fun(c(1,x),beta))
  lines(as.matrix(dat)[,1], preds, ...)
}

#Functions
fun0<-function (x,beta) {plogis(x%*%beta)}  ##logit function
fun1 <- function(x,beta) {exp(x%*%beta)} #Exp Function

#Compute the RMSEs
RMSE_LK_CLK <- function(cubic=F, x1,x2, quadratic=F,pp, fun=fun){
  
  testx <- sort(seq(x1, x2, length.out=1000) )
  
  truth <- pp(testx)
  
  if(cubic){
    test <- data.frame(x=testx,xx=testx*testx,xxx=testx*testx*testx)
  }
  else if(quadratic){
    test <- data.frame(x=testx, xx=testx*testx)
  }
  else{
    test <- data.frame(x=testx)
  } 

    CLK <- CLK[2:ncol(CLK)]
    LK <- LK[2:ncol(LK)]
   
    CLK <- as.matrix(CLK)
    LK <- as.matrix(LK)
    
    for(i in 1:1000){
    #RMSEs
    preds1 <- apply(as.matrix(test), 1, function(x) fun(c(1,x),LK[i,]))
    preds2 <- apply(as.matrix(test), 1, function(x) fun(c(1,x),CLK[i,]))
    
    RMSE <- matrix(0,nrow = 1000,ncol = 2)
    RMSE[i,1]<-sqrt(mean(preds1-truth)^2)
    RMSE[i,2]<-sqrt(mean(preds2-truth)^2)
  }
  #return(colMeans(RMSE))
  write.csv(colMeans(RMSE), "RMSE_Category_1.csv")
}
#Read the saved beta file for the relavent function
CLK <- read.csv("CLK_filename.csv")
LK <- read.csv("LK_filename.csv")

#Category 1 Small size 
#Linear
RMSE_LK_CLK(cubic=F, quadratic=F,x1=-1,x2=1,pp=logistic_3, fun=fun0)

#Quadratic
RMSE_LK_CLK(cubic=F, quadratic=T,x1=-1,x2=1, pp=Qlogistic_3, fun=fun0)

#Cubic
RMSE_LK_CLK(cubic=T, quadratic=F, x1=-1,x2=1, pp=Clogistic_4, fun=fun0)
