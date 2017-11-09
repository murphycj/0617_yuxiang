
basedir = "$baseDir"
params.seqpy = '/Users/charlesmurphy/Desktop/Research/lib/seqpy/bin'
params.count_file = '/Users/charlesmurphy/Desktop/Research/0617_yuxiang/results/RNAseq_1117/HTSeqCount/HTSeq.gene-symbols.counts.csv'

 deseq2_comparisons = Channel
  .from(
    [
      'P-S-Met_vs_P-R-Met',
      'P-S-Met',
      'P-R-Met',
      ['P-S-Met_1','P-S-Met_2','P-S-Met_3'],
      ['P-R-Met_1','P-R-Met_2','P-R-Met_3']
    ],
    [
      'P-S-HCY_vs_P-R-HCY',
      'P-S-HCY',
      'P-R-HCY',
      ['P-S-HCY_1','P-S-HCY_2','P-S-HCY_3'],
      ['P-R-HCY_1','P-R-HCY_2','P-R-HCY_3']
    ],
    [
      'T-S-Met_vs_T-R-Met',
      'T-S-Met',
      'T-R-Met',
      ['T-S-Met_1','T-S-Met_2','T-S-Met_3'],
      ['T-R-Met_1','T-R-Met_2','T-R-Met_3']
    ],
    [
      'T-S-HCY_vs_T-R-HCY',
      'T-S-HCY',
      'T-R-HCY',
      ['T-S-HCY_1','T-S-HCY_2','T-S-HCY_3'],
      ['T-R-HCY_1','T-R-HCY_2','T-R-HCY_3']
    ],
  )


process run_deseq2{

  storeDir "${baseDir}/${name}"

  input:
    set name, p1, p2, group1, group2 from deseq2_comparisons

  output:
    file("${name}") into deseq2_out_pdf
    set name, file("${name}/*csv") into deseq2_out

  """
  Rscript ${params.seqpy}/deseq2.R \
    -counts ${params.count_file} \
    -G1 ${group1.join(",")} \
    -G2 ${group2.join(",")} \
    -minCount 2 \
    -minRowCount 2 \
    -phenotype ${p1},${p2} \
    -outDir ${name}

  """
}

process run_rnk {

  storeDir "${baseDir}/${name}"

  input:
    set name, file(deseq2_out_file) from deseq2_out

  output:
    set file("*.log2FC.rnk") into rnk_out_log2FC
    set file("*statistic.rnk") into rnk_out_statistic
    set file("*logPvalSignFC.rnk") into rnk_out_logPvalSignFC

  """
  python ~/Desktop/Research/lib/seqpy/bin/deseq2rnk.py \
    --infile ${deseq2_out_file} \
    --ranking log2FC \
    --out ${name}.log2FC.rnk
  """
}
