#This code is for plotting
#logistic paramteric functions fitted with logistic functions
######################################################################################

#Define the interested range of x
x <- seq(-1, 1, by=0.01) # predictor for category 1 & 3
x1 <- seq(0, 1, by=0.01) # predictor for category 2

## plot the response curves: 
## from the solymos and Lele (2015) - LK Method

pdf("fig-logit.pdf", width=8, height=12)
op <- par(mfrow=c(3,1), tcl=-0.5, las=1, cex.main=1.25, cex.lab=1, cex.axis=0.75, mai=c(0.5,0.6,0.5,0.6))

#Category 1
plot(0, type="n", ylim=c(0,1), xlim=c(-xlim,xlim),
     xlab="x", ylab="Probability", main="RSPF & LC Satisfied")
for (i in 1:2)
lines(x, plogis(0.606 - 3.64*x), col="black", lwd=3)
lines(x, plogis(0.606 - 3.64*x + 1.26*x^2), col="red", lwd=3)
lines(x, plogis(-0.284 -2.298*x- 0.232*x^2 + 7.269*x^3), col="green", lwd=3)

#Category 2
plot(0, type="n", ylim=c(0,1), xlim=c(0,xlim),
     xlab="x", ylab="Probability", main="Only LC Satisfied (RSPF Not satisfied)")
lines(x1, 0.99*exp(-(4*x1-2)^2), col="black", lwd=3)
lines(x1, 8.3*plogis(-4+2*x1), col="red", lwd=3)
lines(x1, exp(-0.5 + 0.5*x1), col="green", lwd=3)

#Category 3
plot(0, type="n", ylim=c(0,1), xlim=c(-xlim,xlim),
     xlab="x", ylab="Probability", main="Only RSPF Satisfied (LC not Satisfied)")
lines(x, plogis(0.5855 + 1.064*x), col="black", lwd=3)
lines(x, plogis(0.5855 + 1.064*x- 0.218*x^2), col="red", lwd=3)
lines(x, plogis(0.5855 + 1.064*x- 0.218*x^2- 1.81*x^3), col="green", lwd=3)

par(op)
dev.off()
###################################################################
###################################################################

