#This code shows how to plot the saved beta values of the functions to 
#recreate the graphs in the paper

#Function to add lines 
makeline <- function(beta, dat,fun, ...)
{
  preds <- apply(as.matrix(dat), 1, function(x) fun(c(1,x),beta))
  lines(as.matrix(dat)[,1], preds, ...)
}

fun0<-function (x,beta) {plogis(x%*%beta)}  ##logit function
fun1 <- function(x,beta) {exp(x%*%beta)} #Exp Function

plot_LK_CLK <- function(cubic=F, quadratic=F,pp,xl=T,yl=T, fun=fun){
  
  testx <- sort(seq(-1,1, length.out=1000) )
  
  truth <- pp(testx)
  
  plot(testx, truth, col="black", type="l", xlab=if (xl) "x" else "", ylab="P(y=1|x)", ylim=c(0,1), lwd=2, xaxt='n', main = "Linear Logsitic")
  axis(1, labels=xl)
  
  if (cubic){
    test <- data.frame(x=testx,xx=testx*testx,xxx=testx*testx*testx)
  }
  else if (quadratic){
    test <- data.frame(x=testx, xx=testx*testx)
  }
  else{
    test <- data.frame(x=testx)
  }
  
 for(i in 1:1000){
    
    if (cubic){
      yCLK <- plogis(CLK[i,2] + CLK[i,3]*testx + CLK[i,4]*(testx*testx) + CLK[i,5]*(testx*testx*testx))
      yLK <- plogis(LK[i,2] + LK[i,3]*testx + LK[i,4]*(testx*testx)+ LK[i,5]*(testx*testx*testx))
    }
    else if (quadratic){
      yCLK <- plogis(CLK[i,2] + CLK[i,3]*testx + CLK[i,4]*(testx*testx))
      yLK <- plogis(LK[i,2] + LK[i,3]*testx + LK[i,4]*(testx*testx))
    }               
    else{
      yCLK <- plogis(CLK[i,2] + CLK[i,3]*testx)
      yLK <- plogis(LK[i,2] + LK[i,3]*testx )
      
    }
    lines(testx,yLK,col="#d9ef8b",lty=2)
    lines(testx,yCLK,col = "#fee0b6" ,lty=2)
 }
  LK_N <- as.data.frame(LK,nrow = 1000,ncol = ncol(LK))[2:ncol(LK)]
  CLK_N <- as.data.frame(CLK,nrow = 1000,ncol = ncol(CLK))[2:ncol(LK)]
  
  #To show the robustness of the change in local knowledge
  CLK_L <- as.data.frame(CLK_L,nrow = 1,ncol = ncol(CLK_L))[2:ncol(CLK_L)]
  CLK_U <- as.data.frame(CLK_U,nrow = 1,ncol = ncol(CLK_U))[2:ncol(CLK_U)]
  
  makeline(colMeans(LK_N), test,fun, col="#1b7837", lwd = 2)#LK Method
  makeline(colMeans(CLK_N), test,fun, col="#b35806", lwd = 2)#CLK Method
  
  makeline(colMeans(CLK_L), test,fun, col="brown", lwd = 1,lty=2)#CLK Method LL
  makeline(colMeans(CLK_U), test,fun, col="brown", lwd = 1,lty=2)#CLK Method UL
  
  lines(testx,truth,col="black") 
}

#########################################################
#Figure B.1
pdf("fig_Combined Graph of large functions_Only LK.pdf", width=8, height=6)
op <- par(mfrow=c(1,3), tcl=-0.5, las=1, cex.main=1.25, cex.lab=1, cex.axis=0.75, mai=c(0.5,0.6,0.5,0.6))

setwd("H:/PHD/Machine Learning/Recreating Lele's Functions Using LK Method/Graphs for Paper/PBL Removed/Saved Beta Excel Files/logit large/")

LK <- read.csv("LK_filename_large.csv")
#Linear
plot_LK_CLK(cubic=F, quadratic=F,pp=logistic_1, fun=fun0)
#Quadratic
plot_LK_CLK(cubic=F, quadratic=T,pp=Qlogistic_1, fun=fun0)
#Cubic
plot_LK_CLK(cubic=T, quadratic=F,pp=Clogistic_1, fun=fun0)
par(op)
dev.off()

##########################################################

#Figure C.2
#Cloglog 
pdf("fig_Combined Graph of Clog Functions_RSPF Only with local knowledge.pdf", width=8, height=6)
op <- par(mfrow=c(1,3), tcl=-0.5, las=1, cex.main=1.25, cex.lab=1, cex.axis=0.75, mai=c(0.5,0.6,0.5,0.6))

setwd("C:/Users/s3679640/OneDrive - RMIT University/PHD/Machine Learning/Recreating Lele's Functions Using LK Method/Graphs for Paper/PBL Removed/Saved Beta Excel Files/")

CLK <- read.csv("CLK_filename with_c=0.7.csv")
LK <- read.csv("LK_filename with_c=0.7.csv")
#Linear
plot_LK_CLK(cubic=F, quadratic=F,pp= cloglog_5, fun=fun0)
#Quadratic
plot_LK_CLK(cubic=F, quadratic=T,pp=Qcloglog_2, fun=fun0)
#Cubic
plot_LK_CLK(cubic=T, quadratic=F,pp=Ccloglog_new, fun=fun0)
par(op)
dev.off()

################################################################
#Figure C.3
pdf("fig_Combined Graph of Clog Functions_RSPF and LC.pdf", width=8, height=6)
op <- par(mfrow=c(1,3), tcl=-0.5, las=1, cex.main=1.25, cex.lab=1, cex.axis=0.75, mai=c(0.5,0.6,0.5,0.6))

CLK <- read.csv("CLK_filename.csv")
LK <- read.csv("LK_filename.csv")
#Linear
plot_LK_CLK(cubic=F, quadratic=F,pp= cloglog_1, fun=fun0)
#Quadratic
plot_LK_CLK(cubic=F, quadratic=T,pp=Qcloglog_1, fun=fun0)
#Cubic
plot_LK_CLK(cubic=T, quadratic=F,pp=Ccloglog_4, fun=fun0)
par(op)
dev.off()

######################################

#Figure 3 
pdf("fig_Combined Graph of functions that satisfy both LC & RSPF.pdf", width=8, height=6)
op <- par(mfrow=c(1,3), tcl=-0.5, las=1, cex.main=1.25, cex.lab=1, cex.axis=0.75, mai=c(0.5,0.6,0.5,0.6))

CLK <- read.csv("CLK_filename.csv")
LK <- read.csv("LK_filename.csv")

#Linear
plot_LK_CLK(cubic=F, quadratic=F,pp=logistic_3, fun=fun0)

#Quadratic
plot_LK_CLK(cubic=F, quadratic=T,pp=Qlogistic_3, fun=fun0)

#Cubic
plot_LK_CLK(cubic=T, quadratic=F,pp=Clogistic_4, fun=fun0)
par(op)
dev.off()

###########################################################
#Figure 4 top row
pdf("fig_Combined Graph of functions that satisfy only LC .pdf", width=8, height=6)
op <- par(mfrow=c(1,3), tcl=-0.5, las=1, cex.main=1.25, cex.lab=1, cex.axis=0.75, mai=c(0.5,0.6,0.5,0.6))

CLK <- read.csv("CLK_filename.csv")
LK <- read.csv("LK_filename.csv")
#Linear
plot_LK_CLK(cubic=F, quadratic=T,pp=gaussian, fun=fun1)

#Quadratic
plot_LK_CLK(cubic=F, quadratic=F,pp=semilogistic, fun=fun1)

#Cubic
plot_LK_CLK(cubic=F, quadratic=F,pp=exponen, fun=fun1)
par(op)
dev.off()

############################################################
#Figure 4 bottom row
#category 3 with local knowledge and +- 10% of the c.

pdf("fig_Combined Graph of functions Category 3 with robusness.pdf", width=8, height=6)
op <- par(mfrow=c(1,3), tcl=-0.5, las=1, cex.main=1.25, cex.lab=1, cex.axis=0.75, mai=c(0.5,0.6,0.5,0.6))

CLK <- read.csv("CLK_Filename_with_c=0.83.csv")
LK <- read.csv("LK_Filename_with_c=0.83.csv")

#Qith upper and lower bound
CLK_L <- read.csv("CLK_Filename_With_LB_c=0.73.csv")
CLK_U <- read.csv("CLK_Filename_With_UB_c=0.93.csv")

#Linear
plot_LK_CLK(cubic=F, quadratic=F,pp=logistic_1, fun=fun0)
#Quadratic
plot_LK_CLK(cubic=F, quadratic=T,pp=Qlogistic_1, fun=fun0)
#Cubic
plot_LK_CLK(cubic=T, quadratic=F,pp=Clogistic_1, fun=fun0)
par(op)
dev.off()

#######################
#Figure 7 top row
pdf("fig_Combined Graph of Hastie functions with changing C.pdf", width=8, height=6)
op <- par(mfrow=c(1,2), tcl=-0.5, las=1, cex.main=1.25, cex.lab=1, cex.axis=0.75, mai=c(1,1,1,0.6))

#Scaled function
CLK <- read.csv("CLK_Hastie_half_c=0.4.csv")
LK <- read.csv("LK_Hastie_half_c=0.4.csv")
plot_LK_CLK(cubic=F, quadratic=F,pp=hastie_half, fun=fun0)

#Full function
CLK <- read.csv("CLK_Hastie_full_c=0.8.csv")
LK <- read.csv("LK_Hastie_full_c=0.8.csv")
plot_LK_CLK(cubic=F, quadratic=F,pp=hastie, fun=fun0)

par(op)
dev.off()

##################################################

#Figure 7 bottom row
#Hastie scaled function
pdf("fig_Combined Graph of scaled Hastie function.pdf", width=8, height=6)
op <- par(mfrow=c(1,3), tcl=-0.5, las=1, cex.main=1.25, cex.lab=1, cex.axis=0.75, mai=c(0.5,0.6,0.5,0.6))

CLK <- read.csv("filename.csv")
LK <- read.csv("LK_filename.csv")
#With np=50
plot_LK_CLK(cubic=F, quadratic=F,pp=hastie_half, fun=fun0)

#With np=500
plot_LK_CLK(cubic=F, quadratic=F,pp=hastie_half, fun=fun0)

#With np=5000
plot_LK_CLK(cubic=F, quadratic=F,pp=hastie_half, fun=fun0)

par(op)
dev.off()
###########################################
