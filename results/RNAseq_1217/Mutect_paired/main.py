import vcf

bam_dir = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/0617_yuxiang/results/RNAseq_1217/STAR/'
p_bam = bam_dir + 'P/P.gatk.bam'
p_bambai = bam_dir + 'P/P.gatk.bam.bai'

fout = open('comparison.csv','w')
for sample in ['R1-HCY1','R1-HCY2','R1-HCY3','R1-MET1','R1-MET2','R1-MET3','R2-HCY1','R2-HCY2','R2-HCY3','R2-MET1','R2-MET2','R2-MET3','R1','R2']:
    fout.write(
        sample + ',P,' +
        bam_dir + sample + '/' + sample + '.gatk.bam,' +
        bam_dir + sample + '/' + sample + '.gatk.bam.bai,' +
        p_bam + ',' +
        p_bambai + '\n'
    )

fout.close()
