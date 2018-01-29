nohup ~/chm2059/lib/samtools-1.3.1/samtools merge \
  ./P/P.gatk.bam \
  ./P-HCY1/P-HCY1.gatk.bam \
  ./P-HCY2/P-HCY2.gatk.bam \
  ./P-HCY3/P-HCY3.gatk.bam \
  ./P-MET1/P-MET1.gatk.bam \
  ./P-MET2/P-MET2.gatk.bam \
  ./P-MET3/P-MET3.gatk.bam &

nohup ~/chm2059/lib/samtools-1.3.1/samtools merge \
  ./R1/R1.gatk.bam \
  ./R1-HCY1/R1-HCY1.gatk.bam \
  ./R1-HCY2/R1-HCY2.gatk.bam \
  ./R1-HCY3/R1-HCY3.gatk.bam \
  ./R1-MET1/R1-MET1.gatk.bam \
  ./R1-MET2/R1-MET2.gatk.bam \
  ./R1-MET3/R1-MET3.gatk.bam &

nohup ~/chm2059/lib/samtools-1.3.1/samtools merge \
  ./R2/R2.gatk.bam \
  ./R2-HCY1/R2-HCY1.gatk.bam \
  ./R2-HCY2/R2-HCY2.gatk.bam \
  ./R2-HCY3/R2-HCY3.gatk.bam \
  ./R2-MET1/R2-MET1.gatk.bam \
  ./R2-MET2/R2-MET2.gatk.bam \
  ./R2-MET3/R2-MET3.gatk.bam &

nohup ~/chm2059/lib/samtools-1.3.1/samtools index ./P/P.gatk.bam &
nohup ~/chm2059/lib/samtools-1.3.1/samtools index ./R1/R1.gatk.bam &
nohup ~/chm2059/lib/samtools-1.3.1/samtools index ./R2/R2.gatk.bam &
