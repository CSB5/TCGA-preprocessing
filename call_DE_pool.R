#!/mnt/software/bin/Rscript-3.1.0 --slave

#argv[1]: (normalized read) counts file
#argv[2]: number of normal samples
#argv[3]: file extension for results (default in oncoIMPACT: result)
#argv[4]: Software to run: foldChange or DESeq (!!!! not supported)

library(DESeq2)
argv <- commandArgs(TRUE)
runID <- sub(".dat", "", argv[1])

counts <- read.table(argv[1], header=T)

switch(argv[4], 
	foldChange={
		
		#Construct the dds data
		numNormal <- as.numeric(argv[2])
		numSamples <- numNormal + 1
		colData <- as.data.frame(matrix(ncol=1, nrow=numSamples))
		colnames(colData) <- c("condition")
		rownames(colData) <- colnames(counts)
		colData[,1] <- c(rep("untreated",numNormal), rep("treated",1))
		temp <- apply(counts, 2, as.integer)
		rownames(temp) <- rownames(counts)
		temp <- na.omit(temp)
		dds <- DESeqDataSetFromMatrix(countData = temp, colData = colData, design = ~ condition)

		#Log transform
		rld <- rlogTransformation( dds )
		
		#Get the median control read count
		if(numNormal != 1) median_control=as.matrix(apply(assay(rld)[,c(1:numNormal)], 1, FUN = median)) else median_control = as.matrix(assay(rld)[,c(1:numNormal)])
				
		#Compute the log fold change
		res <- data.frame( assay(rld), 
		                   medianControl = median_control[,1], 
				   avgLogExpr = ( median_control[,1] + assay(rld)[,numSamples] ) / 2, 
				   rLogFC = assay(rld)[,numSamples] - median_control[,1] )
				     
  		write.table(res, file=paste(runID, argv[3], sep="."), sep="\t", row.names=T, col.names=T, quote=FALSE)		   

	},
	DESeq={
		# DESeq - typical run
		
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
		normFactors <- matrix(1, ncol=ncol(dds), nrow=nrow(dds))
		normalizationFactors(dds) <- normFactors
		dds <- estimateDispersions(dds)
		dds <- nbinomWaldTest(dds)

		
		# Perform DE test
		res <- results(dds)
		res <- res[order(res$padj),]
		
		write.table(res, file=paste(runID, argv[3], sep="."), sep="\t", row.names=T, col.names=T, quote=FALSE)
		
	}
)