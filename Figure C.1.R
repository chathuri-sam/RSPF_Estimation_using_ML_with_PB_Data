#This code is for plotting
#cloglog paramteric functions fitted with logistic functions
######################################################################################

#Define the interested range of x
x <- seq(-1, 1, by=0.01) # predictor for category 1 & 3

#Part 2 - Cloglog functions for Mis-specified models
pdf("fig-Clog.pdf", width=8, height=12)
xlim <- 1
op <- par(mfrow=c(2,1), tcl=-0.5, las=1, cex.main=1.25, cex.lab=1, cex.axis=0.75, mai=c(1,1.1,1,1.1))

#Category 1
plot(0, type="n", ylim=c(0,1), xlim=c(-xlim,xlim),
     xlab="x", ylab="Probability", main="RSPF & LC Satisfied")
for (i in 1:2)
  lines(x, binomial("cloglog")$linkinv(0.5855 + 1.064*x), col="black", lwd=3)
lines(x, binomial("cloglog")$linkinv(0.5855 + 1.064*x- 0.218*x^2), col="red", lwd=3)
lines(x, binomial("cloglog")$linkinv(-0.284 -2.298*x- 0.232*x^2 + 7.269*x^3), col="green", lwd=3)

#Category 2
plot(0, type="n", ylim=c(0,1), xlim=c(-xlim,xlim),
     xlab="x", ylab="Probability", main="Only RSPF Satisfied (LC not Satisfied)")
lines(x, binomial("cloglog")$linkinv(-0.886 - 1.16*x), col="black", lwd=3)
lines(x, binomial("cloglog")$linkinv(0.37 + 1.56*x - 1.5*x^2), col="red", lwd=3)
lines(x, binomial("cloglog")$linkinv(0.1 - 0.064*x - 0.85*x^2 - 0.81*x^3), col="green", lwd=3)

par(op)
dev.off()