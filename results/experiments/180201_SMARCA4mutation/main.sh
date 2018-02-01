~/chm2059/lib/vcftools_0.1.13/bin/vcftools --vcf ../../RNAseq_1117/Mutect_paired/mutect.vcf --chr 19 --from-bp 11132513 --to-bp 11132513 --recode

for i in 'T-R-Met_2' 'P-R-HCY_3' 'T-R-Met_3' 'P-R-Met_3' 'T-R-HCY_1' 'T-R-Met_1' 'P-R-Met_2' 'T-R-HCY_2' 'P-R-HCY_2' 'T-R-HCY_3' 'P-R-HCY_1' 'P-R-Met_1'
do
  samtools mpileup \
    -f /athena/elementolab/scratch/chm2059/from_dat02/chm2059/data/refdata/hg19_gatk/Homo_sapiens_assembly19.fasta \
    -r 19:11132513-11132513 \
    ../../RNAseq_1117/STAR/$i/$i.gatk.bam > $i.mpileup
done;

python ~/chm2059/lib/seqpy/bin/aggregate_mpileup.py \
  --vcf out.recode.vcf \
  --pileups P-R-HCY_2.mpileup P-R-Met_1.mpileup P-R-Met_3.mpileup T-R-HCY_2.mpileup T-R-Met_1.mpileup T-R-Met_3.mpileup P-R-HCY_1.mpileup P-R-HCY_3.mpileup P-R-Met_2.mpileup T-R-HCY_1.mpileup T-R-HCY_3.mpileup T-R-Met_2.mpileup \
  --samples P-R-HCY_2 P-R-Met_1 P-R-Met_3 T-R-HCY_2 T-R-Met_1 T-R-Met_3 P-R-HCY_1 P-R-HCY_3 P-R-Met_2 T-R-HCY_1 T-R-HCY_3 T-R-Met_2
