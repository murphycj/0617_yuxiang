library(calibrate)
library(gplots)
library(genefilter)
library(dendextend)
library(xlsx)
library("beeswarm")

perform_pcas <- function(data,sample_data) {

	re <- prcomp(t(data))

	pcs=re$sdev**2
	pcs = pcs/sum(pcs) * 100

	pdf("PC-percent-variance.pdf")
	par(mar=c(5.1, 4.1, 4.1, 4.1), xpd=TRUE,cex=1.2)
	barplot(
		pcs,
		names.arg=unlist(lapply(1:length(pcs),function(x) return(paste("PC",x,sep="")))),
		las=2,
		ylab="Percentage of variance (%)",
		ylim=c(0,100)
	)
	dev.off()

	group_labels = list()
	group_labels[["parental"]] <- 1
	group_labels[["YZ_R"]] <- 2
	group_labels[["QL_R"]] <- 3

	samples <- as.character(sample_data$line)
	samples[samples=="parental"]<-1
	samples[samples=="YZ_R"]<-2
	samples[samples=="QL_R"]<-3
	samples <- as.numeric(samples)

	colors <- rainbow(max(samples),start=0.1,end=0.9)[samples]

	legend_colors = rainbow(max(samples),start=0.1,end=0.9)

	#save loadings
	write.csv(re$rotation,"PCA-loadings.csv")

	#pca
	for (i in 1:NUMBER_PCS) {
		for (j in 1: NUMBER_PCS) {
			if (i == j) {
				next
			}
			pdf(paste("./PCA/PC",i,"PC",j,".pdf",sep=""))
			par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
			plot(
				re$x[,i],
				re$x[,j],
				pch=16,
				col=colors,
				ylab=paste("PC",j,sep=""),
				xlab=paste("PC",i,sep=""),
				main="PCA on all samples",
				ylim=c(1.2*min(re$x[,j]),1.2*max(re$x[,j])),
				xlim=c(1.2*min(re$x[,i]),1.2*max(re$x[,i]))
			)
			legend(x="topright",names(group_labels),cex=0.85,inset=c(-0.35,0),pch=16,col=legend_colors,title="Sample type")
			textxy(re$x[,i],re$x[,j],names(re$x[,i]))
			dev.off()
		}
	}
}


NUMBER_PCS = 6
data <- read.csv("..//HTSeqCount/HTSeq.geneSymbols.counts.vst.csv",row.names=1,check.names=F)

sample_data <- read.xlsx("../../../data/samples.xlsx","samples",check.names=F)

data <- data[,as.character(sample_data$name)]

pdf("mean-sd.pdf")
smoothScatter(apply(data,1,mean),apply(data,1,sd))
dev.off()

perform_pcas(data=data,sample_data=sample_data)

#cluster on correlation distance

diss <- 1 - cor(data,use="complete.obs",method="spearman")
diss.2 <- as.dist(diss)
dend <- as.dendrogram(hclust(diss.2,method="average"))

#original labels
png("./Cluster/cluster-all-samples-correlation-distance.png",width=1200,height=600)
plot(dend,main="Hierarchical clustering on all genes",ylab="Correlation distance")
dev.off()

#clsuter on euclidean distance

diss <- dist(t(data))
dend <- as.dendrogram(hclust(diss,method="average"))

png("./Cluster/cluster-all-samples-euclidean.png",width=1200,height=600)
plot(dend,main="Hierarchical clustering on all genes",ylab="Euclidean distance")
dev.off()
