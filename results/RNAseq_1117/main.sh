Rscript ~/Desktop/Research/lib/seqpy/bin/biomart_convert_gene_identifiers.R \
  -data genes-fpkm.csv -d csv -intype ensembl_gene_id -outtype hgnc_symbol -dataset hsapiens_gene_ensembl -out genes-symbols-fpkm.csv
Rscript ~/Desktop/Research/lib/seqpy/bin/biomart_convert_gene_identifiers.R \
  -data HTSeq.gene.counts.csv -d csv -intype ensembl_gene_id -outtype hgnc_symbol -dataset hsapiens_gene_ensembl -out HTSeq.gene-symbols.counts.csv
