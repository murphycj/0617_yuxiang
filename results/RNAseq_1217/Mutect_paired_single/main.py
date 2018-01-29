import vcf

bam_dir = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/0617_yuxiang/results/RNAseq_1217/STAR/'

fout = open('comparison.csv','w')
for sample in ['R1-HCY1','R1-HCY2','R1-HCY3','R1-MET1','R1-MET2','R1-MET3','R2-HCY1','R2-HCY2','R2-HCY3','R2-MET1','R2-MET2','R2-MET3']:
    tmp = sample.split('-')[1]

    reference = 'P-' + tmp
    fout.write(
        sample + ',' + reference + ',' +
        bam_dir + sample + '/' + sample + '.gatk.bam,' +
        bam_dir + sample + '/' + sample + '.gatk.bam.bai,' +
        bam_dir + reference + '/' + reference + '.gatk.bam,' +
        bam_dir + reference + '/' + reference + '.gatk.bam.bai\n'
    )

fout.close()
