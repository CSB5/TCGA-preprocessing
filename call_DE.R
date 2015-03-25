#!/mnt/software/bin/Rscript-3.1.0 --slave

#argv[1]: counts file
#argv[2]: number of normal samples
#argv[3]: file extension for results (default in oncoIMPACT: result)
#argv[4]: Software to run: edgeR or DESeq


argv <- commandArgs(TRUE)
counts <- read.table(argv[1], header=T)
runID <- sub(".dat", "", argv[1])

switch(argv[4], 
	edgeR={
		# edgeR
		library(edgeR)
		conds <- c(rep("N", argv[2]), rep("T", 1))
		counts <- apply(counts, 2, as.integer)
		counts <- na.omit(counts);
		libSize <- apply(counts, 2, sum)
		d <- DGEList(counts=counts, group=conds, lib.size=libSize)
		d <- calcNormFactors(d, method="upperquartile")
		d <- estimateCommonDisp(d)
		de.com <- exactTest(d)
		options(digits=4)
		detags.com <- rownames(topTags(de.com)$table)
		all <- topTags(de.com, n=nrow(counts))$table
		good <- sum(all$FDR < 0.05)
		goodList <- topTags(de.com, n=good)
		all <- all[order(all$FDR),]
		
		write.table(all, file=paste(runID, argv[3], sep="."), sep="\t", row.names=T, col.names=T, quote=FALSE)
	},
	DESeq={
		# DESeq - typical run
		library(DESeq2)
		numSamples <- as.numeric(argv[2]) + 1	# argv[2]: number of normals
		colData <- as.data.frame(matrix(ncol=1, nrow=numSamples))
		colnames(colData) <- c("condition")
		rownames(colData) <- colnames(counts)
		colData[,1] <- c(rep("untreated",argv[2]), rep("treated",1))
		temp <- apply(counts, 2, as.integer)
		rownames(temp) <- rownames(counts)
		temp <- na.omit(temp)
		dds <- DESeqDataSetFromMatrix(countData = temp, colData = colData, design = ~ condition)
		colData(dds)$condition <- factor(colData(dds)$condition,levels=c("untreated","treated"))
		dds <- DESeq(dds)

		# generate normalized counts file
		write.table(counts(dds, normalized=TRUE), file=paste(runID, paste(argv[4], "normalized_counts.txt", sep="."), sep="."), sep="\t", row.names=T, col.names=T, quote=FALSE)
		
		# Perform DE test
		res <- results(dds)
		res <- res[order(res$padj),]
		
		write.table(res, file=paste(runID, argv[3], sep="."), sep="\t", row.names=T, col.names=T, quote=FALSE)
		
	}
)