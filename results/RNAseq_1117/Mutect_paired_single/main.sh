nohup ~/chm2059/lib/gatk-4.0.0.0/gatk \
  AddOrReplaceReadGroups \
  --INPUT ../STAR/R1/R1.gatk.bam \
  --OUTPUT R1.gatk.bam \
  --RGLB NA \
  --RGPL illumina \
  --RGPU NA \
  --RGSM R1-test &

nohup ~/chm2059/lib/gatk-4.0.0.0/gatk \
  AddOrReplaceReadGroups \
  --INPUT ../STAR/P/P.gatk.bam \
  --OUTPUT P.gatk.bam \
  --RGLB NA \
  --RGPL illumina \
  --RGPU NA \
  --RGSM P &


nohup ~/chm2059/lib/gatk-4.0.0.0/gatk \
  AddOrReplaceReadGroups \
  --INPUT ../STAR/R2/R2.gatk.bam \
  --OUTPUT R2.gatk.bam \
  --RGLB NA \
  --RGPL illumina \
  --RGPU NA \
  --RGSM R2 &

nohup ~/chm2059/lib/gatk-4.0.0.0/gatk \
  Mutect2 \
  --reference /athena/elementolab/scratch/chm2059/from_dat02/chm2059/data/refdata/hg19_gatk/Homo_sapiens_assembly19.fasta \
  --input R1.gatk.bam \
  --input P.gatk.bam \
  --tumor-sample R1 \
  --normal-sample P \
  --create-output-variant-index \
  --min-base-quality-score 15 \
  --output R1.vcf \
  --native-pair-hmm-threads 12 \
  --minimum-mapping-quality 20 &

nohup ~/chm2059/lib/gatk-4.0.0.0/gatk \
  Mutect2 \
  --reference /athena/elementolab/scratch/chm2059/from_dat02/chm2059/data/refdata/hg19_gatk/Homo_sapiens_assembly19.fasta \
  --input R2.gatk.bam \
  --input P.gatk.bam \
  --tumor-sample R2 \
  --normal-sample P \
  --create-output-variant-index \
  --min-base-quality-score 15 \
  --output R2.vcf \
  --native-pair-hmm-threads 12 \
  --minimum-mapping-quality 20 &

~/chm2059/lib/gatk-4.0.0.0/gatk FilterMutectCalls \
  --variant R1.vcf \
  --output tmp.vcf \
  --normal-artifact-lod -1 \
  --tumor-lod 8

awk -F '\t' '{if($0 ~ /\#/) print; else if($7 == "PASS") print}' tmp.vcf > R1.tmp.vcf

~/chm2059/lib/gatk-4.0.0.0/gatk FilterMutectCalls \
  --variant R2.vcf \
  --output tmp.vcf \
  --normal-artifact-lod -2 \
  --base-quality-score-threshold 10 \
  --use-new-qual-calculator \
  --tumor-lod 12

awk -F '\t' '{if($0 ~ /\#/) print; else if($7 == "PASS") print}' tmp.vcf > R2.tmp.vcf

~/chm2059/lib/bcftools-1.3.1/bcftools view -Oz -s R1 -o R1.filtered.vcf.gz -I R1.tmp.vcf
~/chm2059/lib/bcftools-1.3.1/bcftools view -Oz -s R2 -o R2.filtered.vcf.gz -I R2.tmp.vcf

rm *tmp*

/athena/elementolab/scratch/chm2059/from_dat02/chm2059/lib/tabix-0.2.6/tabix R1.filtered.vcf.gz
/athena/elementolab/scratch/chm2059/from_dat02/chm2059/lib/tabix-0.2.6/tabix R2.filtered.vcf.gz

export PATH=~/chm2059/lib/tabix-0.2.6/:$PATH
/athena/elementolab/scratch/chm2059/from_dat02/chm2059/lib/bcftools-1.3.1/bcftools merge -m none R1.filtered.vcf.gz R2.filtered.vcf.gz > mutect.vcf

/athena/elementolab/scratch/chm2059/from_dat02/chm2059/lib/jre1.8.0_25/bin/java -Xmx16g -jar \
  /athena/elementolab/scratch/chm2059/from_dat02/chm2059/lib/snpEff/SnpSift.jar annotate \
  /athena/elementolab/scratch/chm2059/from_dat02/chm2059/data/refdata/dbsnp/human_9606_b149_GRCh37p13/VCF/common_all_20161121.vcf.gz \
  mutect.vcf \
  | /athena/elementolab/scratch/chm2059/from_dat02/chm2059/lib/jre1.8.0_25/bin/java -Xmx16g -jar \
  /athena/elementolab/scratch/chm2059/from_dat02/chm2059/lib/snpEff/snpEff.jar eff \
  -v GRCh37.75 -canon \
  | /athena/elementolab/scratch/chm2059/from_dat02/chm2059/lib/jre1.8.0_25/bin/java -Xmx16g -jar \
  /athena/elementolab/scratch/chm2059/from_dat02/chm2059/lib/snpEff/SnpSift.jar filter \
  "! exists ID" > mutect.annotated.vcf

python /athena/elementolab/scratch/chm2059/from_dat02/chm2059/lib/seqpy/bin/vcf_to_table.py \
  --vcf mutect.annotated.vcf \
  --out mutect.annotated.csv \
  --mutect

python /athena/elementolab/scratch/chm2059/from_dat02/chm2059/lib/seqpy/bin/vcf_to_table.py \
  --vcf mutect.annotated.vcf \
  --out nodbsnp.mutect.csv \
  --everything \
  --mutect


~/chm2059/lib/gatk-4.0.0.0/gatk \
  AddOrReplaceReadGroups \
  --INPUT R1.gatk.bam \
  --OUTPUT R1-test.gatk.bam \
  --RGLB NA \
  --RGPL illumina \
  --RGPU NA \
  --RGSM R1-test
