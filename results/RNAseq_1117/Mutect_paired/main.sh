nohup java -jar ~/chm2059/lib/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar \
  -T MuTect2 \
  -R /athena/elementolab/scratch/chm2059/from_dat02/chm2059/data/refdata/hg19_gatk/Homo_sapiens_assembly19.fasta \
  -I:tumor ../STAR/T-R/T-R.gatk.bam \
  -I:normal ../STAR/T-S/T-S.gatk.bam \
  --min_base_quality_score 15 \
  --max_alt_alleles_in_normal_count 3 \
  -o T-R.vcf \
  -nct 12 &

nohup java -jar ~/chm2059/lib/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar \
  -T MuTect2 \
  -R /athena/elementolab/scratch/chm2059/from_dat02/chm2059/data/refdata/hg19_gatk/Homo_sapiens_assembly19.fasta \
  -I:tumor ../STAR/P-R/P-R.gatk.bam \
  -I:normal ../STAR/P-S/P-S.gatk.bam \
  --min_base_quality_score 15 \
  --max_alt_alleles_in_normal_count 3 \
  -o P-R.vcf \
  -nct 12 &
