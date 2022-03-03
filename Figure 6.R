#Plots the RMSEs of all the categories together as given in Figure 6

pdf("Bar_Chart_RMSEs with Small sample.pdf",  width=8, height=6)
op <- par(mfrow=c(2,2), tcl=-0.5, las=1, cex.main=1.25, cex.lab=1, cex.axis=0.75, mai=c(0.5,0.8,0.5,0.6))

#Revised LC and RSPF
dat1 <- c(0.00002340243,	0.000018508020,
          0.000006888289,	0.000022195480,
          0.000001208790,	0.000005828634)

# Generate the bar chart for Figure 1
dat <- t(matrix(dat1, ncol=3))
dat <- dat[c(1,2,3),]

barx <- barplot(t(dat[,c(1,2)]), beside=T, ylim=c(0,0.0008),args.legend = list(bty = "n", x = "top", ncol = 1), legend.text=c("LK","CLK"), names.arg=c("Linear Logit", "Quad logit", "Cubic logit"), ylab="RMS error",main="Category 1")
########################################
#For the LC Only Functions (SMALL SAMPLE)
dat2 <- c(0.0000451253800,	0.0000132746600,
         0.0001132478,	0.0000628082600,	
         0.00005400429, 0.000032009970)


# Generate the bar chart for Figure 1
dat <- t(matrix(dat2, ncol=3))
dat <- dat[c(1,2,3),]

barx <- barplot(t(dat[,c(1,2)]), beside=T, ylim=c(0,0.0008),args.legend = list(bty = "n", x = "top", ncol = 1), legend.text=c("LK","CLK"), names.arg=c("Linear Logit", "Quad logit", "Cubic logit"), ylab="RMS error",main="Category 1 (Small Sample)")

#######################
#For the LC Only Functions 
dat3 <- c(0.00043333240, 0.0000034553860,
         0.0004360376,	0.0000028621270,	
         0.00061809720,	0.0000201624800)


# Generate the bar chart for Figure 1
dat <- t(matrix(dat3, ncol=3))
dat <- dat[c(1,2,3),]

barx <- barplot(t(dat[,c(1,2)]), beside=T, ylim=c(0,0.0008),args.legend = list(bty = "n", x = "top", ncol = 1), legend.text=c("LK","CLK"), names.arg=c("SemiLogit", "Gaussian","Exponential" ), ylab="RMS error",main="Category 2")

##################################
#For Only RSPF Function with (10%) different local knowledge
dat4 <- c(0.0001012236000,	0.0000289518400,
          0.0005232283000,	0.0000293827300,
          0.0005214721000,	0.0000032844240)

# Generate the bar chart for Figure 1
dat <- t(matrix(dat4, ncol=3))
dat <- dat[c(1,2,3),]

barx <- barplot(t(dat[,c(1,2)]), beside=T, ylim=c(0,0.0008),args.legend = list(bty = "n", x = "top", ncol = 1), legend.text=c("LK","CLK"), names.arg=c("Linear logit" , "Quad logit", "Cubic logit"), ylab="RMS error",main="Category 3")
dev.off()

