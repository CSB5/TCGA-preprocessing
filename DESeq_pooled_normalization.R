#!/mnt/software/bin/Rscript-3.1.0 --slave

#argv[1]: counts file
#argv[2]: number of normal samples

argv <- commandArgs(TRUE)
counts <- read.table(argv[1], header=T, check.names=F)
runID <- sub(".dat", "", argv[1])

rowSum <- apply(counts, 1, sum)
counts <- counts[rowSum > 0,] + 1


library(DESeq2)
numSamples <- dim(counts)[2]
colData <- as.data.frame(matrix(ncol=1, nrow=numSamples))
colnames(colData) <- c("condition")
rownames(colData) <- colnames(counts)
colData[,1] <- c(rep("untreated",argv[2]), rep("treated",(numSamples - as.integer(argv[2]))))
temp <- apply(counts, 2, as.integer)
rownames(temp) <- rownames(counts)
temp <- na.omit(temp)
dds <- DESeqDataSetFromMatrix(countData = temp, colData = colData, design = ~ condition)
colData(dds)$condition <- factor(colData(dds)$condition,levels=c("untreated","treated"))
dds <- DESeq(dds)
			
# generate normalized counts file
write.table(counts(dds, normalized=TRUE), file=paste(runID, "DESeq.normalized_counts.txt", sep="."), sep="\t", row.names=T, col.names=T, quote=FALSE)
save(dds, file=paste(runID, "DESeq.normalized_counts.rda", sep="."))
