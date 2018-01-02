date=180102

ln /Users/charlesmurphy/Desktop/Research/0617_yuxiang/results/RNAseq_1217/Cufflinks/genes-symbols-fpkm.csv $date.genes-symbols-fpkm.csv
cp /Users/charlesmurphy/Desktop/Research/0617_yuxiang/results/RNAseq_1217/STAR/STAR.QC.csv $date.STAR.QC.csv
cp /Users/charlesmurphy/Desktop/Research/0617_yuxiang/results/RNAseq_1217/UnsupervisedAnalyses/FPKM/PCA/all/PC1PC2.pdf $date.PC1PC2.pdf
cp /Users/charlesmurphy/Desktop/Research/0617_yuxiang/results/RNAseq_1217/UnsupervisedAnalyses/FPKM/Cluster/cluster-all-samples-euclidean.png $date.cluster-all-samples-euclidean.png

for i in 'PHcy_vs_R1Hcy' 'PHcy_vs_R2Hcy' 'PHcy_vs_RHcy' 'PMet_vs_PHcy' 'PMet_vs_R1Met' 'PMet_vs_R2Met' 'PMet_vs_RMet' 'P_vs_R' 'R1Met_vs_R1Hcy' 'R2Met_vs_R2Hcy'
do
  cp /Users/charlesmurphy/Desktop/Research/0617_yuxiang/results/RNAseq_1217/DESeq2/$i/$i\_results.csv $date.$i\_results.csv
done;
