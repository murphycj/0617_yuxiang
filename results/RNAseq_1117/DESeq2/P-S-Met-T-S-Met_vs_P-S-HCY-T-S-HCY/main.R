
library("DESeq2")
library("RColorBrewer")
library("gplots")
library("genefilter")

group1 <- c('P-S-Met_1','P-S-Met_2','P-S-Met_3','T-S-Met_1','T-S-Met_2','T-S-Met_3')
group2 <- c('T-S-HCY_1','T-S-HCY_2','T-S-HCY_3','P-S-HCY_1','P-S-HCY_2','P-S-HCY_3')
phenotype <- c("P-S-Met-T-S-Met","P-S-HCY-T-S-HCY")
comparison <- "P-S-Met-T-S-Met_vs_P-S-HCY-T-S-HCY"

sink(paste("samples_compared.txt",sep=""))
cat("Group1\n")
cat(paste(group1,collapse=", "))
cat("\n")
cat("Group2\n")
cat(paste(group2,collapse=", "))
cat("\n")
sink()

data <- read.csv("/Users/charlesmurphy/Desktop/Research/0617_yuxiang/results/RNAseq_1117/HTSeqCount/HTSeq.gene-symbols.counts.csv", header=T, row.names=1, check.names=FALSE)
data <- data[,c(group1, group2)]
data <- data[rowSums(data)>=2,]

coldata <- data.frame(
  mainFactor=factor(c(
    rep(phenotype[1],length(group1)),
    rep(phenotype[2],length(group2))
    ),
    levels=phenotype
  ),
  cellline=factor(c(rep("P",3),rep("T",3),rep("T",3),rep("P",3)),levels=c("T","P")),
  row.names=colnames(data)
)

countTable <- DESeqDataSetFromMatrix(
  countData=data,
  colData=coldata,
  design=~cellline+mainFactor
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
