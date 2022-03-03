#Plots the RMSEs of the two categories of Cloglog functions
dat <- c(0.000083208300,	0.000056705010,
         0.000117332100,	0.000011240310,
         0.000156505,	0.000145597)

# Generate the bar chart for Figure 1
dat <- t(matrix(dat, ncol=3))
dat <- dat[c(1,2,3),]

pdf("Bar_Chart_RMSEs_1000iterations_Cloglog Only RSPF Satisfied.pdf",  width=8, height=6)
barx <- barplot(t(dat[,c(1,2)]), beside=T, ylim=c(0,0.0009),args.legend = list(bty = "n", x = "top", ncol = 1), legend.text=c("LK","CLK"), names.arg=c("Linear Clog", "Quadratic Clog", "Cubic Clog"), ylab="RMS error")
dev.off()

#For the RSPF and LC Functions PBL Removed
dat <- c(0.000078747,	0.000018820,
         0.000677795,	0.000002644,	
         0.00001298867,	0.00000631770)

# Generate the bar chart for Figure 1
dat <- t(matrix(dat, ncol=3))
dat <- dat[c(1,2,3),]

pdf("Bar_Chart_RMSEs_1000iterations_Cloglog RSPF and LC Satisfied Selected.pdf",  width=8, height=6)
barx <- barplot(t(dat[,c(1,2)]), beside=T, ylim=c(0,0.0008),args.legend = list(bty = "n", x = "top", ncol = 1), legend.text=c("LK","CLK"), names.arg=c("Linear Clog", "Quadratic Clog", "Cubic Clog"), ylab="RMS error")
dev.off()
