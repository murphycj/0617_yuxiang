# git:

export date="063017"

cp ../../results/RNAseq/ClusteringVST/PCA/*
cp ../../results/RNAseq/ClusteringVST/Cluster/*
cp ../../results/RNAseq/Cufflinks/filtered/genesSymbols-fpkm.csv .

find *txt | parallel mv {} $date.{}
