library(circlize)
library(xlsx)

de_genes <- read.csv(
  "../../RNAseq/DESeq2/parental_vs_resistant_ensembl/parental_vs_resistant_results.csv",
  header=T,
  check.names=F,
  row.names=1,
  colClasses=c("character","numeric","numeric","numeric","numeric","numeric","numeric")
)
de_genes <- de_genes[(!is.na(de_genes$padj)) & (de_genes$padj<=0.01) & (abs(de_genes$log2FoldChange)>=1),]

gene_pos <- read.csv("ensembl75.csv",row.names=1,header=F,colClasses=c("character","character","numeric","numeric","character"))
gene_pos <- gene_pos[row.names(de_genes),]

colors <- rep("red",length(de_genes$log2FoldChange))
colors[de_genes$log2FoldChange<0]<-"blue"

results <- read.xlsx("./parental_vs_resistant/parental_vs_resistant_results_up.xlsx","Sheet1",colClasses=c(rep("numeric",7),"character"))
results <- results[results$num_differentially_expressed_genes>3,]
results$chrom <- unlist(lapply(results$chrom,function(x) return(paste("chr",x,sep=""))))

de_data <- data.frame(
  chr=gene_pos[,1],
  start=gene_pos[,2],
  end=gene_pos[,3],
  value1=de_genes$log2FoldChange,
  color=rep("red",length(de_genes$log2FoldChange)),
  stringsAsFactors=FALSE
)
de_data[de_data$value1<0,"color"]<-"blue"
de_data$chr <- unlist(lapply(de_data$chr,function(x) return(paste("chr",x,sep=""))))
de_data <- de_data[de_data$chr!="chrMT",]


pdf("parental_vs_resistant_results_circos.pdf")
circos.par("gap.degree"=c(rep(1,23),7),"start.degree"=90)
circos.initializeWithIdeogram(
  species="hg19"
)
circos.genomicTrackPlotRegion(
  de_data[,c("chr","start","end","value1")]
)
circos.yaxis(labels.cex = 0.6,sector.index="chr1")
for (chrom in unique(de_data$chr)) {
  tmp <- de_data[de_data$chr==chrom,]
  circos.genomicPoints(region=tmp[,c("start","end")],value=data.frame(value=tmp$value1),col=tmp$color,sector.index=chrom,pch=".",cex=2.5)
}
circos.genomicTrackPlotRegion(
  data.frame(
    chr=results$chrom,
    start=results$start,
    end=results$end,
    stringsAsFactors=FALSE
  ),
  ylim=c(0,65)
)
circos.yaxis(labels.cex = 0.6,sector.index="chr1")
for (chrom in unique(results$chrom)) {
  tmp <- results[results$chrom==chrom,]
  circos.genomicLines(
    region=tmp[,c("start","end")],
    value=data.frame(value=-log10(tmp$padj)),
    type = "segment",
    col="red",
    sector.index=chrom,
    lwd = 4
  )
}
dev.off()


circos.clear()
pdf("parental_vs_resistant_results_circos.chr1.pdf",width=20,height=20)
circos.par("gap.degree"=c(330),"start.degree"=c(90))
circos.initializeWithIdeogram(
  species="hg19",
  chromosome.index="chr1"
)
circos.genomicTrackPlotRegion(
  de_data[,c("chr","start","end","value1")],
  ylim=c(-4,4)
)
circos.yaxis(labels.cex = 0.6,sector.index="chr1")
chrom="chr1"
tmp <- de_data[de_data$chr==chrom,]
tmp$color[6] <- "green"
tmp <- rbind(tmp,tmp[6,])
tmp <- tmp[-c(6),]
circos.genomicPoints(region=tmp[,c("start","end")],value=data.frame(value=tmp$value1),col=tmp$color,sector.index=chrom,pch=".",cex=3.5)

circos.genomicTrackPlotRegion(
  data.frame(
    chr=results$chrom,
    start=results$start,
    end=results$end,
    stringsAsFactors=FALSE
  ),
  ylim=c(0,65)
)
circos.yaxis(labels.cex = 0.6,sector.index="chr1")
chrom="chr1"
tmp <- results[results$chrom==chrom,]
circos.genomicLines(
  region=tmp[,c("start","end")],
  value=data.frame(value=-log10(tmp$padj)),
  type = "segment",
  col="red",
  sector.index=chrom,
  lwd = 4
)
dev.off()
