#!/mnt/software/bin/Rscript-3.1.0 --slave

library(vioplot)
v1 <- read.csv(file="basic_stats.dat", sep="\t",head=TRUE)

postscript(file="diff_supp_0_4.eps",
           paper="special",
           width=10,
           height=10,
           horizontal=FALSE)
vioplot(v1$DYS_0_5, v1$DYS_1, v1$DYS_1_5, v1$DYS_2, v1$DYS_2_5 ,v1$DYS_3, v1$DYS_3_5, v1$DYS_4, ylim = c(0,6500), col = "lightblue", names=c(">=0.5", ">=1", ">=1.5", ">=2", ">=2.5", ">=3", ">=3.5", ">=4"))
title(main = "Absolute Log2 fold change\ncompared to normal gene expression", cex.main = 2)
abline(h=2315, col = "red", lwd = 3, lty=3)
abline(h=1852, col = "blue", lwd = 3, lty=3 )
abline(h=926, col = "green", lwd = 3, lty=3)
abline(h=463, col = "black", lwd = 3, lty=3)

postscript(file="alt_vioplot.eps",
           paper="special",
           width=5,
           height=5,
           horizontal=FALSE)

vioplot(v1$MUT, v1$CNV, (v1$MUT+v1$CNV), col = "indianred", names=c("SNV", "CNV", "SNV+CNV"))

postscript(file="snv_distrib.eps",
           paper="special",
           width=5,
           height=5,
           horizontal=FALSE)
bins=seq(-25, 4000,by=50)
avg = round(mean(v1$MUT), 2)
hist.data = hist(v1$MUT, breaks = bins, plot=F);
plot(hist.data,  col = "indianred", xlab = "nb non synonymous SNVs", main = paste("Avg nb of non synonymous SNVs", avg))

postscript(file="cnv_distrib.eps",
           paper="special",
           width=5,
           height=5,
           horizontal=FALSE)
bins=seq(-25, 2000,by=50)
avg = round(mean(v1$CNV), 2)
hist.data = hist(v1$CNV, breaks = bins, plot=F)
plot(hist.data = "indianred", main = paste("Avg nb CNV gene", avg))

postscript(file="alteration_distrib.eps",
           paper="special",
           width=5,
           height=5,
           horizontal=FALSE)
bins=seq(-25, 5000,by=50)
avg = round(mean(v1$MUT+v1$CNV))
hist.data = hist(v1$MUT+v1$CNV, breaks =bins, plot=F)
#hist.data$counts = (hist.data$counts)/127
plot(hist.data, col = "indianred", xlab = "nb altered genes", main = paste("Avg nb of altered genes", avg))

postscript(file="diff_supp_0.eps",
           paper="special",
           width=10,
           height=10,
           horizontal=FALSE)
vioplot(v1$DYS_0_5, v1$DYS_1, v1$DYS_1_5, v1$DYS_2, v1$DYS_2_5 ,v1$DYS_3, v1$DYS_3_5, v1$DYS_4, v1$DYS_4_5, v1$DYS_5, v1$DYS_5_5, v1$DYS_6, ylim = c(0,6500), col = "lightblue", names=c(">=0.5", ">=1", ">=1.5", ">=2", ">=2.5", ">=3", ">=3.5", ">=4", ">=4.5", ">=5", ">=5.5", ">=6"))
title(main = "Absolute Log2 fold change\ncompared to normal gene expression", cex.main = 2)
abline(h=2315, col = "red", lwd = 3, lty=3)
abline(h=1852, col = "blue", lwd = 3, lty=3 )
abline(h=926, col = "green", lwd = 3, lty=3)
abline(h=463, col = "black", lwd = 3, lty=3)

postscript(file="diff_all.eps",
           paper="special",
           width=5,
           height=5,
           horizontal=FALSE)

vioplot(v1$DYS_0,  v1$DYS_0_5, v1$DYS_1, v1$DYS_1_5, v1$DYS_2, v1$DYS_2_5 ,v1$DYS_3, v1$DYS_3_5, v1$DYS_4, v1$DYS_4_5, v1$DYS_5, v1$DYS_5_5, v1$DYS_6, ylim = c(0,9261), col = "lightblue", names=c(">0", ">=0.5", ">=1", ">=1.5", ">=2", ">=2.5", ">=3", ">=3.5", ">=4", ">=4.5", ">=5", ">=5.5", ">=6"))
title("Absolute Log2 fold change\ncompared to normal gene expression")
