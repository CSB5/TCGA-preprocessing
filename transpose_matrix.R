#!/mnt/software/bin/Rscript-3.1.0 --slave

argv <- commandArgs(TRUE)
data <- read.table(argv[1], header=T, check.names=F)
write.table(t(data), file=argv[2], quote=F, sep= "\t", row.names=T, col.names=T)
