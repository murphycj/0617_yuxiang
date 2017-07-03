library("DESeq2")
library("RColorBrewer")
library("gplots")
library("genefilter")
library("xlsx")

samples <- read.xlsx("../../../../data/samples.xlsx","samples",check.names=F)
row.names(samples) <- as.character(samples$name)
samples <- samples[samples$line=="parental",]
g1 <- as.character(samples[samples$medium=="Met+Hcy-","name"])
g2 <- as.character(samples[samples$medium=="Met-Hcy+","name"])
samples <- samples[c(g1,g2),]

sink("samples_compared.txt")
cat("Group1\n")
cat(paste(g1,collapse=", "))
cat("\n")
cat("Group2\n")
cat(paste(g2,collapse=", "))
cat("\n")
sink()

data <- read.csv("../../HTSeqCount/HTSeq.geneSymbols.counts.csv", header=T, row.names=1, check.names=FALSE)
data <- data[,c(g1, g2)]
data <- data[rowSums(data)>=2,]

coldata <- data.frame(
  medium=factor(samples$medium_formated,levels=c("MetPlusHcyNeg","MetNegHcyPlus")),
  row.names=samples$name
)

countTable <- DESeqDataSetFromMatrix(
  countData=data,
  colData=coldata,
  design=~medium
)

result <- DESeq(countTable)
res <- results(result)
dd = res[with(res,order(padj)),]

write.csv(dd,paste("Met+Hcy-_vs_Met-Hcy+_results.csv",sep=""))

pdf("Met+Hcy-_vs_Met-Hcy+-pvalue-hist.pdf")
hist(res$pvalue, breaks=20, col="grey")
dev.off()

pdf("Met+Hcy-_vs_Met-Hcy+-padj-hist.pdf")
hist(res$padj, breaks=20, col="grey")
dev.off()

pdf("Met+Hcy-_vs_Met-Hcy+-MA-plot.pdf")
plotMA(res)
dev.off()
