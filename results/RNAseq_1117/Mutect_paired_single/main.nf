#!/usr/bin/env nextflow

basedir = "$baseDir"

params.reference = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/data/refdata/hg19_gatk/Homo_sapiens_assembly19.fasta'
params.reference_dict = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/data/refdata/hg19_gatk/Homo_sapiens_assembly19.dict'
params.reference_fai = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/data/refdata/hg19_gatk/Homo_sapiens_assembly19.fasta.fai'
params.gtf = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/data/refdata/hg19_gatk/Homo_sapiens.GRCh37.85.gtf'
params.dbsnp = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/data/refdata/dbsnp/human_9606_b149_GRCh37p13/VCF/common_all_20161121.vcf.gz'

pairs_mutect = Channel
  .from(file("comparison.csv"))
  .splitCsv(header: false)
  .map{it -> [it[0], it[1], it[2], it[3], it[4], it[5]]}

process mutect {

  storeDir "${baseDir}/${tumor}"
  scratch true
  executor 'sge'
  clusterOptions '-l h_vmem=2G -pe smp 12 -l h_rt=96:00:00 -l athena=true'

  input:
    set tumor, control, tumor_bam, tumor_bam_bai, control_bam, control_bam_bai from pairs_mutect

  output:
    set tumor, file('*vcf') into mutect_paired_out

  """
  rsync -L ${tumor_bam} tumor.bam
  rsync -L ${tumor_bam_bai} tumor.bam.bai
  rsync -L ${control_bam} control.bam
  rsync -L ${control_bam_bai} control.bam.bai

  ~/chm2059/lib/gatk-4.0.0.0/gatk \
    Mutect2 \
    --reference /athena/elementolab/scratch/chm2059/from_dat02/chm2059/data/refdata/hg19_gatk/Homo_sapiens_assembly19.fasta \
    --input tumor.bam \
    --input control.bam \
    --tumor-sample ${tumor} \
    --normal-sample ${control} \
    --create-output-variant-index \
    --min-base-quality-score 15 \
    --output ${tumor}.vcf \
    --native-pair-hmm-threads 12 \
    --minimum-mapping-quality 20
  """
}

process filter_mutect {

  storeDir "${baseDir}/${prefix}"

  input:
    set prefix, vcf from mutect_paired_out

  output:
    file('*.filtered.mutect.vcf.gz') into filter_mutect_out
    file('*.filtered.mutect.vcf.gz.tbi') into filter_mutect_out_tbi

  """
  ~/chm2059/lib/gatk-4.0.0.0/gatk FilterMutectCalls \
    --variant ${vcf} \
    --output tmp.vcf \
    --normal-artifact-lod -1 \
    --tumor-lod 8

  python ${params.seqpy}/filter_mutect.py \
    --vcf tmp.vcf \
    --AD 4 \
    --freq 0.05 \
    --min_tumor 12 \
    --min_normal 8 \
    --tlod 8 \
    --tumor ${prefix} \
    --normal P \
    --max_normal_support 3 \
    --out tmp2.vcf

  ${params.vcftools} --vcf tmp2.vcf --remove-indv P2 --recode

  python ${params.seqpy}/modify_QSS_field.py \
    --vcf out.recode.vcf \
    --out tmp3.vcf

  ${params.vcftools} --vcf tmp3.vcf --indv ${prefix} --recode

  mv out.recode.vcf ${prefix}.filtered.mutect.vcf

  ${params.bgzip} ${prefix}.filtered.mutect.vcf
  ${params.tabix} ${prefix}.filtered.mutect.vcf.gz
  """
}

vcfs_mutect = filter_mutect_out.toList()
vcfs_mutect_tbi = filter_mutect_out_tbi.toList()

process combine_mutect {

  storeDir "${baseDir}"

  input:
    file vcfs_mutect
    file vcfs_mutect_tbi

  output:
    file('mutect.vcf') into combine_mutect_out2
    file('nodbsnp.mutect.vcf') into combine_mutect_out

  """
  export PATH=~/chm2059/lib/tabix-0.2.6/:$PATH

  ${params.bcftools} merge -m none ${vcfs_mutect} > mutect.vcf

  ${params.java} -Xmx16g -jar \
    ${params.snpsift} annotate \
    ${params.dbsnp} \
    mutect.vcf \
    | ${params.java} -Xmx16g -jar \
    ${params.snpeff} eff \
    -v GRCm38.86 -canon \
    | ${params.java} -Xmx16g -jar \
    ${params.snpsift} filter \
    "! exists ID" > tmp2.vcf

  ${params.vcftools} \
    --vcf tmp2.vcf \
    --bed ${params.bed_file} \
    --recode --recode-INFO-all

  mv out.recode.vcf nodbsnp.mutect.vcf
  """
}

process vcf_to_table_mutect {

  storeDir "${baseDir}"

  input:
    file(combine_mutect_out_nodbsnp) from combine_mutect_out

  output:
    set file('nodbsnp.mutect.csv'), file('nodbsnp.mutect.allEffects.csv') into vcf_to_table_out

  """
  python ${params.seqpy}/vcf_to_table.py \
    --vcf ${combine_mutect_out_nodbsnp} \
    --out nodbsnp.mutect.csv \
    --mutect

  python ${params.seqpy}/vcf_to_table.py \
    --vcf ${combine_mutect_out_nodbsnp} \
    --everything \
    --out nodbsnp.mutect.allEffects.csv \
    --mutect
  """
}
