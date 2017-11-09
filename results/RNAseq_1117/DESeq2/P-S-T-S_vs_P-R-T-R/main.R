
library("DESeq2")
library("RColorBrewer")
library("gplots")
library("genefilter")


group1 <- c('T-S-Met_1','T-S-Met_2','T-S-Met_3','T-S-HCY_1','T-S-HCY_2','T-S-HCY_3','P-S-Met_1','P-S-Met_2','P-S-Met_3','P-S-HCY_1','P-S-HCY_2','P-S-HCY_3')
group2 <- c('T-R-Met_1','T-R-Met_2','T-R-Met_3','T-R-HCY_1','T-R-HCY_2','T-R-HCY_3','P-R-Met_1','P-R-Met_2','P-R-Met_3','P-R-HCY_1','P-R-HCY_2','P-R-HCY_3')
phenotype <- c("P-S-T-S","P-R-T-R")
comparison <- "P-S-T-S_vs_P-R-T-R"

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
  media=factor(c("Met","Met","Met","HCY","HCY","HCY","Met","Met","Met","HCY","HCY","HCY","Met","Met","Met","HCY","HCY","HCY","Met","Met","Met","HCY","HCY","HCY"),levels=c("Met","HCY")),
  cellline=factor(c(rep("T",6),rep("P",6),rep("T",6),rep("P",6)),levels=c("T","P")),
  row.names=colnames(data)
)

countTable <- DESeqDataSetFromMatrix(
  countData=data,
  colData=coldata,
  design=~cellline+media+mainFactor
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
