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
    set tumor, control, file('*vcf') into mutect_out

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

  storeDir "${baseDir}/${tumor}"

  input:
    set tumor, control, vcf from mutect_out

  output:
    file('*.filtered.mutect.vcf.gz') into filter_mutect_out
    file('*.filtered.mutect.vcf.gz.tbi') into filter_mutect_out_tbi
    file('*.filteredWithNormal.mutect.vcf') into filter_mutect_out2

  """
  ~/chm2059/lib/gatk-4.0.0.0/gatk FilterMutectCalls \
    --variant ${vcf} \
    --output tmp.vcf \
    --normal-artifact-lod 0 \
    --max-events-in-region 4 \
    --tumor-lod 8

  ${params.bcftools} norm -m- -N tmp.vcf > tmp2.vcf

  python ${params.seqpy}/filter_mutect_v4.py \
    --vcf tmp2.vcf \
    --AD 6 \
    --freq 0.1 \
    --min_tumor 12 \
    --min_normal 8 \
    --tlod 8 \
    --nlod 8 \
    --tumor ${tumor} \
    --normal ${control} \
    --max_normal_af 0.05 \
    --out tmp3.vcf

  awk -F '\\t' '{if(\$0 ~ /\\#/) print; else if(\$7 == "PASS") print}' tmp3.vcf > ${tumor}.filteredWithNormal.mutect.vcf

  ${params.bcftools} view -Ov -s ${tumor} ${tumor}.filteredWithNormal.mutect.vcf > ${tumor}.filtered.mutect.vcf

  ${params.bgzip} ${tumor}.filtered.mutect.vcf
  ${params.tabix} ${tumor}.filtered.mutect.vcf.gz
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
    file('mutect.annotated.vcf') into combine_mutect_out

  """
  export PATH=~/chm2059/lib/tabix-0.2.6/:$PATH

  ${params.bcftools} merge -m none ${vcfs_mutect} > mutect.vcf

  ${params.java} -Xmx16g -jar \
    ${params.snpsift} annotate \
    ${params.dbsnp} \
    mutect.vcf \
    | ${params.java} -Xmx16g -jar \
    ${params.snpeff} eff \
    -v GRCh37.75 -canon \
    | ${params.java} -Xmx16g -jar \
    ${params.snpsift} filter \
    "! exists ID" > mutect.annotated.vcf
  """
}

process vcf_to_table_mutect {

  storeDir "${baseDir}"

  input:
    file(combine_mutect_out_nodbsnp) from combine_mutect_out

  output:
    set file('mutect.annotated.csv'), file('mutect.annotated.allEffects.csv') into vcf_to_table_out

  """
  python ${params.seqpy}/vcf_to_table.py \
    --vcf ${combine_mutect_out_nodbsnp} \
    --out mutect.annotated.csv \
    --mutect

  python ${params.seqpy}/vcf_to_table.py \
    --vcf ${combine_mutect_out_nodbsnp} \
    --everything \
    --out mutect.annotated.allEffects.csv \
    --mutect
  """
}
