#!/mnt/software/bin/Rscript-3.1.0 --slave
# argv[1] = counts matrix
# argv[2] = excluded_samples.txt
# argv[3] = output filename

argv <- commandArgs(TRUE)
data <- read.table(argv[1], header=T, check.names=F)
excluded <- read.table(argv[2], header=F, check.names=F)
selected <- setdiff(colnames(data), unlist(excluded))
filtered <- subset(data, select=selected)
write.table(filtered, file=argv[3], quote=F, sep= "\t", row.names=T, col.names=T)