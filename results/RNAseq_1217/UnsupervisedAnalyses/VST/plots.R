library(calibrate)
library(gplots)
library(genefilter)
library(dendextend)
library(xlsx)
library("beeswarm")
library(RColorBrewer)

perform_pcas <- function(data,sample_data,outdir) {

	dir.create(paste("PCA/",outdir,sep=""),showWarnings=F)

	re <- prcomp(t(data))

	pcs=re$sdev**2
	pcs = pcs/sum(pcs) * 100

	pdf(paste("./PCA/",outdir,"/PC-percent-variance.pdf",sep=""))
	par(mar=c(5.1, 4.1, 4.1, 4.1), xpd=TRUE,cex=1.2)
	barplot(
		pcs,
		names.arg=unlist(lapply(1:length(pcs),function(x) return(paste("PC",x,sep="")))),
		las=2,
		ylab="Percentage of variance (%)",
		ylim=c(0,100)
	)
	dev.off()

	groups <- unique(as.character(sample_data$group))

	colors_ref <- brewer.pal(length(groups),"Spectral")
	colors <- as.character(sample_data$group)
	for (i in 1:length(groups)) {
		colors[colors==groups[i]] <- colors_ref[i]
	}

	#save loadings
	write.csv(re$rotation,paste("./PCA/",outdir,"/PCA-loadings.csv",sep=""))

	#pca
	for (i in 1:NUMBER_PCS) {
		for (j in 1: NUMBER_PCS) {
			if (i == j) {
				next
			}
			pdf(paste("./PCA/",outdir,"/PC",i,"PC",j,".pdf",sep=""))
			par(mar=c(5.1, 4.1, 4.1, 10.1), xpd=TRUE)
			plot(
				re$x[,i],
				re$x[,j],
				pch=16,
				col=colors,
				ylab=paste("PC",j," (",round(pcs[j],1),"%)",sep=""),
				xlab=paste("PC",i," (",round(pcs[i],1),"%)",sep=""),
				main="",
				ylim=c(1.7*min(re$x[,j]),1.7*max(re$x[,j])),
				xlim=c(1.7*min(re$x[,i]),1.7*max(re$x[,i]))
			)
			legend(x="topright",groups,cex=0.85,inset=c(-0.45,0),pch=16,col=colors_ref,title="Sample type")
			textxy(re$x[,i],re$x[,j],names(re$x[,i]))
			dev.off()
		}
	}
}


NUMBER_PCS = 6
data <- read.csv("/Users/charlesmurphy/Desktop/Research/0617_yuxiang/results/RNAseq_1217/HTSeqCount/HTSeq.gene-symbols.vst.csv",row.names=1,check.names=F)

samples <- read.xlsx("/Users/charlesmurphy/Desktop/Research/0617_yuxiang/data/samples.xlsx","Sheet1",strinsAsFactors=F)
samples <- samples[as.character(samples$sequence_date)=="Yuxiang5316_2017_12_22",]
samples$group <- paste(samples$cell_line,samples$media_conditions,sep="-")

data <- data[,as.character(samples$sample_ID)]

pdf("mean-sd.pdf")
smoothScatter(apply(data,1,mean),apply(data,1,sd))
dev.off()

perform_pcas(data=data,sample_data=samples,outdir="all")

data.sub <- data[order(apply(data,1,var),decreasing=T),]
perform_pcas(data=data.sub[1:500,],sample_data=samples,outdir="top500Var")
perform_pcas(data=data.sub[1:1000,],sample_data=samples,outdir="top1000Var")
perform_pcas(data=data.sub[1:2000,],sample_data=samples,outdir="top2000Var")

data.sub <- data[order(apply(data,1,mean),decreasing=T),]
perform_pcas(data=data.sub[1:500,],sample_data=samples,outdir="top500Mean")
perform_pcas(data=data.sub[1:1000,],sample_data=samples,outdir="top1000Mean")
perform_pcas(data=data.sub[1:2000,],sample_data=samples,outdir="top2000Mean")
