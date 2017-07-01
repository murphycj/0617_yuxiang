library(DESeq2)
library(xlsx)

samples <- read.xlsx("../../../data/samples.xlsx","samples")

counts <- read.csv("HTSeq.geneSymbols.counts.csv",row.names=1,header=T,check.names=F)
counts <- counts[rowSums(counts)>1,]
counts <- counts[,samples$name]

coldata <- data.frame(
  group=samples$line,
  row.names=colnames(counts)
)

dds <- DESeqDataSetFromMatrix(counts,coldata,~group)

vst <- varianceStabilizingTransformation(dds,blind=T)

write.csv(assay(vst),"HTSeq.geneSymbols.counts.vst.csv",quote=F)
