
library("DESeq2")
library("RColorBrewer")
library("gplots")
library("genefilter")


group1 <- c("P-HCY1","P-HCY2","P-HCY3")
group2 <- c("R2-HCY1","R2-HCY2","R2-HCY3")
phenotype <- c("PHcy","R2Hcy")
comparison <- "PHcy_vs_R2Hcy"
setwd(comparison)

sink(paste("samples_compared.txt",sep=""))
cat("Group1\n")
cat(paste(group1,collapse=", "))
cat("\n")
cat("Group2\n")
cat(paste(group2,collapse=", "))
cat("\n")
sink()

data <- read.csv("/Users/charlesmurphy/Desktop/Research/0617_yuxiang/results/RNAseq_1217/HTSeqCount/HTSeq.gene-symbols.counts.csv", header=T, row.names=1, check.names=FALSE)
data <- data[,c(group1, group2)]
data <- data[rowSums(data)>=2,]

coldata <- data.frame(
  mainFactor=factor(c(
    rep(phenotype[1],length(group1)),
    rep(phenotype[2],length(group2))
    ),
    levels=phenotype
  ),
  row.names=colnames(data)
)

countTable <- DESeqDataSetFromMatrix(
  countData=data,
  colData=coldata,
  design=~mainFactor
)

result <- DESeq(countTable)

res <- results(result)

dd = res[with(res,order(padj)),]

write.csv(
  dd,
  paste(comparison,"_results.csv",sep="")
)

pdf(paste(comparison,"-pvalue-hist.pdf",sep=""))
hist(res$pvalue, breaks=20, col="grey")
dev.off()

pdf(paste(comparison,"-padj-hist.pdf",sep=""))
hist(res$padj, breaks=20, col="grey")
dev.off()

pdf(paste(comparison,"-MA-plot.pdf",sep=""))
plotMA(res)
dev.off()
