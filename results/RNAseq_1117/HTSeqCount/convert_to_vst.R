library(DESeq2)
library(xlsx)

samples <- read.xlsx("/Users/charlesmurphy/Desktop/Research/0617_yuxiang/data/samples.xlsx","Sheet1",strinsAsFactors=F)
samples <- samples[as.character(samples$sequence_date)=="17/11/1",]

counts <- read.csv("HTSeq.gene-symbols.counts.csv",row.names=1,header=T,check.names=F)
counts <- counts[rowSums(counts)>1,]
counts <- counts[,as.character(samples$sample_ID)]

coldata <- data.frame(
  group=samples$cell_line,
  row.names=colnames(counts)
)

dds <- DESeqDataSetFromMatrix(counts,coldata,~group)
vst <- varianceStabilizingTransformation(dds,blind=T)
write.csv(assay(vst),"HTSeq.gene-symbols.vst.csv",quote=F)
