library("DESeq2")
library("RColorBrewer")
library("gplots")
library("genefilter")
library("xlsx")

samples <- read.xlsx("../../../../data/samples.xlsx","samples",check.names=F)
row.names(samples) <- as.character(samples$name)
samples <- samples[samples$medium=="Met+Hcy-",]
g1 <- as.character(samples[samples$line=="parental","name"])
g2 <- as.character(samples[samples$line!="parental","name"])
samples <- samples[c(g1,g2),]

sink("samples_compared.txt")
cat("Group1\n")
cat(paste(g1,collapse=", "))
cat("\n")
cat("Group2\n")
cat(paste(g2,collapse=", "))
cat("\n")
sink()

data <- read.csv("../../HTSeqCount/HTSeq.gene.counts.csv", header=T, row.names=1, check.names=FALSE)
data <- data[,c(g1, g2)]
data <- data[rowSums(data)>=2,]

cell_line <- as.character(samples$line)
cell_line[cell_line=="YZ_R"] <- "resistant"
cell_line[cell_line=="QL_R"] <- "resistant"

coldata <- data.frame(
  line=factor(cell_line,levels=c("parental","resistant")),
  row.names=samples$name
)

countTable <- DESeqDataSetFromMatrix(
  countData=data,
  colData=coldata,
  design=~line
)

result <- DESeq(countTable)
res <- results(result)
dd = res[with(res,order(padj)),]

write.csv(dd,paste("parental_vs_resistant_Met+Hcy-_results.csv",sep=""))

pdf("parental_vs_resistant_Met+Hcy-_pvalue-hist.pdf")
hist(res$pvalue, breaks=20, col="grey")
dev.off()

pdf("parental_vs_resistant_Met+Hcy-_padj-hist.pdf")
hist(res$padj, breaks=20, col="grey")
dev.off()

pdf("parental_vs_resistant_Met+Hcy-_MA-plot.pdf")
plotMA(res)
dev.off()
