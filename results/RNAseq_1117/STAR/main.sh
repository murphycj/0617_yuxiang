#nohup ~/chm2059/lib/samtools-1.3.1/samtools merge \
#  ./P-R/P-R.gatk.bam \
#  ./P-R-HCY_1/P-R-HCY_1.gatk.bam \
#  ./P-R-HCY_2/P-R-HCY_2.gatk.bam \
#  ./P-R-HCY_3/P-R-HCY_3.gatk.bam \
#  ./P-R-Met_1/P-R-Met_1.gatk.bam \
#  ./P-R-Met_2/P-R-Met_2.gatk.bam \
#  ./P-R-Met_3/P-R-Met_3.gatk.bam &

nohup ~/chm2059/lib/samtools-1.3.1/samtools merge \
  ./P-S/P-S.gatk.bam \
  ./P-S-HCY_1/P-S-HCY_1.gatk.bam \
  ./P-S-HCY_2/P-S-HCY_2.gatk.bam \
  ./P-S-HCY_3/P-S-HCY_3.gatk.bam \
  ./P-S-Met_1/P-S-Met_1.gatk.bam \
  ./P-S-Met_2/P-S-Met_2.gatk.bam \
  ./P-S-Met_3/P-S-Met_3.gatk.bam &

#nohup ~/chm2059/lib/samtools-1.3.1/samtools merge \
#  ./T-R/T-R.gatk.bam \
#  ./T-R-HCY_1/T-R-HCY_1.gatk.bam \
#  ./T-R-HCY_2/T-R-HCY_2.gatk.bam \
#  ./T-R-HCY_3/T-R-HCY_3.gatk.bam \
#  ./T-R-Met_1/T-R-Met_1.gatk.bam \
#  ./T-R-Met_2/T-R-Met_2.gatk.bam \
#  ./T-R-Met_3/T-R-Met_3.gatk.bam &

nohup ~/chm2059/lib/samtools-1.3.1/samtools merge \
  ./T-S/T-S.gatk.bam \
  ./T-S-HCY_1/T-S-HCY_1.gatk.bam \
  ./T-S-HCY_2/T-S-HCY_2.gatk.bam \
  ./T-S-HCY_3/T-S-HCY_3.gatk.bam \
  ./T-S-Met_1/T-S-Met_1.gatk.bam \
  ./T-S-Met_2/T-S-Met_2.gatk.bam \
  ./T-S-Met_3/T-S-Met_3.gatk.bam &

#nohup ~/chm2059/lib/samtools-1.3.1/samtools index ./P-R/P-R.gatk.bam &
nohup ~/chm2059/lib/samtools-1.3.1/samtools index ./P-S/P-S.gatk.bam &
#nohup ~/chm2059/lib/samtools-1.3.1/samtools index ./T-R/T-R.gatk.bam &
nohup ~/chm2059/lib/samtools-1.3.1/samtools index ./T-S/T-S.gatk.bam &


nohup ~/chm2059/lib/gatk-4.0.0.0/gatk \
  AddOrReplaceReadGroups \
  --INPUT ./P-S/P-S.gatk.bam \
  --OUTPUT ./P-S/P-S.gatk2.bam \
  --RGLB NA \
  --RGPL illumina \
  --RGPU NA \
  --RGSM P-S &

nohup ~/chm2059/lib/gatk-4.0.0.0/gatk \
  AddOrReplaceReadGroups \
  --INPUT ./T-S/T-S.gatk.bam \
  --OUTPUT ./T-S/T-S.gatk2.bam \
  --RGLB NA \
  --RGPL illumina \
  --RGPU NA \
  --RGSM T-S &

rm ./P-S/P-S.gatk.bam*
rm ./T-S/T-S.gatk.bam*

mv ./P-S/P-S.gatk2.bam ./P-S/P-S.gatk.bam
mv ./T-S/T-S.gatk2.bam ./T-S/T-S.gatk.bam

nohup ~/chm2059/lib/samtools-1.3.1/samtools index ./P-S/P-S.gatk.bam &
nohup ~/chm2059/lib/samtools-1.3.1/samtools index ./T-S/T-S.gatk.bam &
