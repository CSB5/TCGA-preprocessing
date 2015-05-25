#!/mnt/software/bin/Rscript-3.1.0 --slave

M = read.table(file="RNA_SEQ_COUNT_MATRIX.DESeq.normalized_counts.txt", header = TRUE, check.names=F)
M_t = t(log2(M + 1))
dist_M = dist(M_t)

library(ctc)
hc <- hclust(d=dist_M)
write.table(hc2Newick(hc), file="hc.newick", quote=F, col.names=F, row.names=F)

library("gplots")
pdf(file='h_clustering.pdf')
sample_col = c(rep("red", 33), rep("blue", 280))
heatmap.2(as.matrix(dist_M), ColSideColors = sample_col)
dev.off()

d <- as.dist(cor(M, method="spearman"))
hc <- hclust(d=d)
write.table(hc2Newick(hc), file="hc_spearman.newick")

pdf(file='h_spearman_clustering.pdf')
sample_col = c(rep("red", 33), rep("blue", 280))
heatmap.2(as.matrix(d), ColSideColors = sample_col)
dev.off()

