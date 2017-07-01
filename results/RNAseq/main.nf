#!/usr/bin/env nextflow

basedir = "$baseDir"

params.reference = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/data/refdata/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
params.reference_dict = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/data/refdata/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.dict'
params.reference_fai = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/data/refdata/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai'
params.gtf = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/data/refdata/GRCh38/Homo_sapiens.GRCh38.88.chr.gtf'

params.star = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/lib/STAR-2.5.3a/bin/Linux_x86_64_static/STAR'
params.starRef = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/data/refdata/GRCh38/star50'
params.outFilterMismatchNmax = 5
params.star_threads = 6

params.dexseq = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/lib/DEXSeq/python_scripts/dexseq_count.py'
params.fastqscreen = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/lib/fastq_screen_v0.10.0/fastq_screen'
params.fastqscreen_conf = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/lib/fastq_screen_v0.10.0/fastq_screen.conf'
params.java = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/lib/jre1.8.0_25/bin/java'
params.fastqc = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/lib/FastQC/fastqc'
params.seqpy = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/lib/seqpy/bin/'
params.samtools = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/lib/samtools-1.1/samtools'
params.picard = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/lib/picard-tools-1.137/picard.jar'
params.htseq = '~/.local/bin/htseq-count'
params.flagstat = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/lib/seqpy/bin/flagstat.py'
params.cufflinks = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/lib/cufflinks-2.2.1.Linux_x86_64/cufflinks'

fqfiles = Channel.create()
fqfiles2 = Channel.create()
fqfiles3 = Channel.create()
fqfiles4 = Channel.create()

Channel
  .fromFilePairs("/athena/elementolab/scratch/chm2059/from_dat02/chm2059/0617_yuxiang/data/Yuxiang4823_2017_06_23/*_L00*_{R1,R2}_001.fastq.gz", size: 2)
  .map{[it[0].split('_S')[0],it[1][0], it[1][1]]}
  .separate(fqfiles, fqfiles2, fqfiles3, fqfiles4){it->[it,it,it,it]}

process fastqc {

  storeDir "${baseDir}/fastqc/${prefix}"
  maxForks 6

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
  executor 'sge'
  scratch true
  clusterOptions '-l h_vmem=10G -pe smp 6 -l h_rt=96:00:00 -l athena=true'

  storeDir "${baseDir}/STAR/${prefix}"

  input:
    set prefix, file(read1), file(read2) from fqfiles3

  output:
    set prefix, file('*bam'), file('*bam.bai') into star_out
    set prefix, file('*bam'), file('*bam.bai') into star_out2
    set prefix, file('*bam'), file('*bam.bai') into star_out3
    set prefix, file('*bam'), file('*bam.bai') into star_out4
    set file('*Log.progress.out'), file('*Log.out'), file('*Chimeric.out.junction'), file('*Chimeric.out.sam') into star_log
    file('*Log.final.out') into star_log_Final
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

   ${params.samtools} index ${prefix}*.bam
   rm -rf 1pass
   rm -rf star_2pass
 """
}

process samtools_flagstat {

  maxForks 8
  storeDir "${baseDir}/SamtoolsFlagstat/${prefix}"

  input:
    set prefix, file(bam_file), file(bam_index_file) from star_out

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
  executor 'sge'
  clusterOptions '-l h_vmem=3G -pe smp 6 -l h_rt=96:00:00 -l athena=true'

  storeDir "${baseDir}/Cufflinks/${prefix}"

  input:
    set prefix, file(bam_file), file(bam_bai_file) from star_out2
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
    set prefix, file(bam_file), file(bam_index_file) from star_out3

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
    file log_files from star_log_Final.toList()
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
