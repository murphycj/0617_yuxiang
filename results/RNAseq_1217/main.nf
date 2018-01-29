#!/usr/bin/env nextflow

basedir = "$baseDir"

params.reference = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/data/refdata/hg19_gatk/Homo_sapiens_assembly19.fasta'
params.reference_dict = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/data/refdata/hg19_gatk/Homo_sapiens_assembly19.dict'
params.reference_fai = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/data/refdata/hg19_gatk/Homo_sapiens_assembly19.fasta.fai'
params.gtf = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/data/refdata/hg19_gatk/Homo_sapiens.GRCh37.85.gtf'
params.dbsnp = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/data/refdata/dbsnp/human_9606_b149_GRCh37p13/VCF/common_all_20161121.vcf.gz'

params.starRef = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/data/refdata/hg19_gatk/star99'
params.outFilterMismatchNmax = 5
params.star_threads = 5

fqfiles = Channel.create()
fqfiles2 = Channel.create()
fqfiles3 = Channel.create()
fqfiles4 = Channel.create()

Channel
  .fromFilePairs("/athena/elementolab/scratch/chm2059/from_dat02/chm2059/0617_yuxiang/data/Yuxiang5316_2017_12_22/*S*_L00{1,2}_{R1,R2}_001.fastq.gz", size: 4)
  .map{[it[0].split('_S')[0],[it[1][0], it[1][2]], [it[1][1], it[1][3]]]}
  .separate(fqfiles, fqfiles2, fqfiles3, fqfiles4){it->[it,it,it,it]}

process fastqc {

  storeDir "${baseDir}/fastqc/${prefix}"

  input:
    set prefix, file(read1), file(read2) from fqfiles

  output:
    file('*zip') into fastqc_results
    val prefix into samples

  """
  zcat ${read1} ${read2} | ${params.fastqc} stdin --outdir=.
  mv stdin_fastqc.zip ${prefix}.zip
  """
}

process star {

  maxForks 6
  storeDir "${baseDir}/STAR/${prefix}"

  input:
    set prefix, file(read1), file(read2) from fqfiles3

  output:
    set prefix, file('*.gatk.bam'), file('*.gatk.bam.bai') into gatk_out
    set prefix, file('*.gatk.bam'), file('*.gatk.bam.bai') into gatk_out2
    set prefix, file('*.gatk.bam'), file('*.gatk.bam.bai') into gatk_out_insert_size
    set prefix, file('*.gatk.bam'), file('*.gatk.bam.bai') into gatk_out3
    set prefix, file('*.gatk.bam'), file('*.gatk.bam.bai') into gatk_out_mutect
    file('*Log.progress.out') into star_out_progress
    file('*Log.final.out') into star_log
    file('*SJ.out.tab') into star_out_SJ
    val prefix into samples_star_qc


  script:
  read1=read1.toString().replaceAll(/ /,",")
  read2=read2.toString().replaceAll(/ /,",")

  """
  mkdir 1pass
  mkdir star_2pass
  ${params.star} \
    --outSAMtype BAM SortedByCoordinate \
    --runThreadN ${params.star_threads} \
    --outSAMstrandField intronMotif \
    --outFilterMismatchNmax ${params.outFilterMismatchNmax} \
    --genomeDir ${params.starRef} \
    --readFilesCommand zcat \
    --readFilesIn ${read1} ${read2} \
    --outFileNamePrefix ./1pass/tmp.

  ${params.star} \
    --runMode genomeGenerate \
    --genomeDir ./star_2pass \
    --genomeFastaFiles ${params.reference} \
    --sjdbFileChrStartEnd ./1pass/tmp.SJ.out.tab \
    --sjdbOverhang 50 \
    --runThreadN ${params.star_threads}

  ${params.star} \
    --outSAMtype BAM SortedByCoordinate \
    --runThreadN ${params.star_threads} \
    --outSAMstrandField intronMotif \
    --genomeDir ./star_2pass \
    --readFilesCommand zcat \
    --readFilesIn ${read1} ${read2} \
    --outFileNamePrefix ${prefix}. \
    --chimOutType SeparateSAMold \
    --chimSegmentMin 1

   ${params.samtools} index ${prefix}.Aligned.sortedByCoord.out.bam
   rm -rf 1pass
   rm -rf star_2pass

   ${params.java} -Xmx10g -jar \
     ${params.picard} AddOrReplaceReadGroups \
     VALIDATION_STRINGENCY=LENIENT \
     I=${prefix}.Aligned.sortedByCoord.out.bam \
     O=tmp.bam \
     ID=${prefix} \
     LB=NA \
     PL=illumina \
     PU=NA \
     RGSM=${prefix}

   ${params.samtools} index tmp.bam

   ${params.java} -Xmx10g -jar \
     ${params.picard} MarkDuplicates \
     VALIDATION_STRINGENCY=LENIENT \
     I=tmp.bam \
     O=tmp2.bam \
     REMOVE_DUPLICATES=false \
     M=duplicates.bam

   ${params.samtools} index tmp2.bam
   rm tmp.bam*

   ${params.java} -jar \
     ${params.gatk} \
     -T SplitNCigarReads \
     -R ${params.reference} \
     -I tmp2.bam \
     -o tmp.bam \
     -rf ReassignOneMappingQuality \
     -RMQF 255 \
     -RMQT 60 \
     -U ALLOW_N_CIGAR_READS

   ${params.samtools} index tmp.bam
   rm tmp2.bam*

   ${params.java} -Xmx10g -jar \
     ${params.gatk} \
     -T RealignerTargetCreator \
     -R ${params.reference} \
     -I tmp.bam \
     -U ALLOW_N_CIGAR_READS \
     -o forIndelRealigner.intervals \
     -nt 4

   ${params.java} -Xmx10g -jar \
     ${params.gatk} \
     -T IndelRealigner \
     -R ${params.reference} \
     -I tmp.bam \
     -targetIntervals forIndelRealigner.intervals \
     -o tmp2.bam

   ${params.samtools} index tmp2.bam

   ${params.java} -Xmx10g -jar \
     ${params.gatk} \
     -T BaseRecalibrator \
     -R ${params.reference} \
     -I tmp2.bam \
     -o recal_data.table \
     -knownSites ${params.dbsnp} \
     -nct 4

   ${params.java} -Xmx10g -jar \
     ${params.gatk} \
     -T PrintReads \
     -R ${params.reference} \
     -I tmp2.bam \
     --BQSR recal_data.table \
     -o tmp3.bam

   ${params.samtools} index tmp3.bam

   python ${params.seqpy}/filter_bad_cigar.py \
     --infile tmp3.bam \
     --out ${prefix}.gatk.bam

   ${params.samtools} index ${prefix}.gatk.bam
 """
}

process samtools_flagstat {

  maxForks 8
  storeDir "${baseDir}/SamtoolsFlagstat/${prefix}"

  input:
    set prefix, file(bam_file), file(bam_index_file) from gatk_out

  output:
    file('*flagstat') into flagstat_out

  """
  ${params.samtools} flagstat ${bam_file} > ${prefix}.flagstat

  """
}

process combine_flagstat {

  storeDir "${baseDir}/SamtoolsFlagstat"

  input:
    file(flagstat_files) from flagstat_out.toList()
  output:
    file('flagstats.csv') into flagstat_combine_out

  """
  python ${params.seqpy}/aggregate_flagstat.py \
    --flagstat ${flagstat_files} \
    --out flagstats.csv
  """
}

process cufflinks {
  scratch true
  maxForks 6

  storeDir "${baseDir}/Cufflinks/${prefix}"

  input:
    set prefix, file(bam_file), file(bam_bai_file) from gatk_out2
  output:
    set prefix, file('*genes.fpkm_tracking'), file('*isoforms.fpkm_tracking'), file('*skipped.gtf'), file('*transcripts.gtf') into cufflinks_out
    file('*genes.fpkm_tracking') into cufflinks_out_genes
    file('*isoforms.fpkm_tracking') into cufflinks_out_isoforms
    val prefix into samples_cufflinks

  """
  ${params.cufflinks} \
    -q -p 6 \
    -o . \
    -b ${params.reference} \
    -G ${params.gtf} \
    ${bam_file}
  mv genes.fpkm_tracking ${prefix}.genes.fpkm_tracking
  mv isoforms.fpkm_tracking ${prefix}.isoforms.fpkm_tracking
  mv skipped.gtf ${prefix}.skipped.gtf
  mv transcripts.gtf ${prefix}.transcripts.gtf
  """
}

process combine_cufflinks {

  storeDir "${baseDir}/Cufflinks/"

  input:
    file genes from cufflinks_out_genes.toList()
    file isoforms from cufflinks_out_isoforms.toList()
    val all_samples from samples_cufflinks.toList()

  output:
    file('isoforms-fpkm.csv') into combine_cufflinks_isoforms
    file('genes-fpkm.csv') into combine_cufflinks_genes

  """
  python /athena/elementolab/scratch/chm2059/from_dat02/chm2059/lib/seqpy/bin/aggregate_fpkm.py \
    --files ${genes} \
    --samples ${all_samples.join(" ")} \
    --duplicates highest \
    --out genes-fpkm.csv

  python /athena/elementolab/scratch/chm2059/from_dat02/chm2059/lib/seqpy/bin/aggregate_fpkm.py \
    --files ${isoforms} \
    --samples ${all_samples.join(" ")} \
    --duplicates highest \
    --out isoforms-fpkm.csv

  """
}

process htseq_reads {

  storeDir "${baseDir}/HTSeqCount/${prefix}"

  input:
    set prefix, file(bam_file), file(bam_index_file) from gatk_out3

  output:
    file('*count') into htseq_reads_out
    val prefix into samples_combine_htseq

  """
  python ${params.htseq} \
    -s no \
    -f bam \
    ${bam_file} \
    ${params.gtf} > ${prefix}.count
  """
}

process combine_htseq {

  storeDir "${baseDir}/HTSeqCount/"

  input:
    file results from htseq_reads_out.toList()
    val all_samples from samples_combine_htseq.toList()

  output:
    file('HTSeq.gene.counts.csv') into combine_htseq_out

  """
  python ${params.seqpy}/aggregate_htseqcount.py \
    --files ${results} \
    --samples ${all_samples.join(" ")} \
    --out HTSeq.gene.counts.csv
  """
}

process combine_star_qc {

  storeDir "${baseDir}/STAR/"

  input:
    file log_files from star_log.toList()
    val all_samples from samples_star_qc.toList()

  output:
    file('STAR.QC.csv') into star_qc_out

  """
  python ${params.seqpy}/aggregate_star_QC.py \
    --files ${log_files} \
    --samples ${all_samples.join(" ")} \
    --out STAR.QC.csv
  """
}

process mutect {

  storeDir "${baseDir}/Mutect/${prefix}"
  scratch true
  executor 'sge'
  clusterOptions '-l h_vmem=2G -pe smp 12 -l h_rt=96:00:00 -l athena=true'

  input:
    set prefix, file(bam_file), file(bam_index_file) from gatk_out_mutect

  output:
    set prefix, file('*vcf') into mutect_out

  """
  rsync -L ${bam_file} tmp.bam
  rsync -L ${bam_index_file} tmp.bam.bai

  ~/chm2059/lib/gatk-4.0.0.0/gatk \
    Mutect2 \
    --reference /athena/elementolab/scratch/chm2059/from_dat02/chm2059/data/refdata/hg19_gatk/Homo_sapiens_assembly19.fasta \
    --input tmp.bam \
    --tumor-sample ${prefix} \
    --create-output-variant-index \
    --min-base-quality-score 15 \
    --output ${prefix}.vcf \
    --native-pair-hmm-threads 12 \
    --minimum-mapping-quality 20
  """
}

process filter_mutect {

  storeDir "${baseDir}/Mutect/${prefix}"

  input:
    set prefix, vcf from mutect_out

  output:
    file('*.filtered.vcf.gz') into filter_mutect_out
    file('*.filtered.vcf.gz.tbi') into filter_mutect_out_tbi

  """

  ~/chm2059/lib/gatk-4.0.0.0/gatk FilterMutectCalls \
    --variant ${vcf} \
    --output tmp.vcf \
    --base-quality-score-threshold 10 \
    --tumor-lod 10 \
    --max-events-in-region 4

  python ${params.seqpy}/filter_mutect.py \
    --vcf tmp.vcf \
    --tumor ${prefix} \
    --nonpaired \
    --AD 6 \
    --freq 0.05 \
    --min_tumor 12 \
    --tlod 10 \
    --out tmp2.vcf

  awk -F '\t' '{if(\$0 ~ /\\#/) print; else if(\$7 == "germline_risk") print}' tmp2.vcf > ${prefix}.filtered.vcf

  ${params.bgzip} ${prefix}.filtered.vcf
  ${params.tabix} ${prefix}.filtered.vcf.gz
  """
}

vcfs_mutect = filter_mutect_out.toList()
vcfs_mutect_tbi = filter_mutect_out_tbi.toList()

process combine_mutect {

  storeDir "${baseDir}/Mutect"

  input:
    file vcfs_mutect
    file vcfs_mutect_tbi

  output:
    set file('mutect.vcf'), file('mutect.annotated.vcf') into combine_mutect_out2

  """
  export PATH=/athena/elementolab/scratch/chm2059/from_dat02/chm2059/lib/tabix-0.2.6/:$PATH

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
