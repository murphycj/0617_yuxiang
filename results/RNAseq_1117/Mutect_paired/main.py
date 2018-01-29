import vcf

bam_dir = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/0617_yuxiang/results/RNAseq_1117/STAR/'
p_bam = bam_dir + 'P-S/P-S.gatk.bam'
p_bambai = bam_dir + 'P-S/P-S.gatk.bam.bai'

fout = open('comparison.csv','w')
for sample in ['T-R-HCY_1','T-R-HCY_2','T-R-HCY_3','T-R-Met_1','T-R-Met_2','T-R-Met_3','P-R-HCY_1','P-R-HCY_2','P-R-HCY_3','P-R-Met_1','P-R-Met_2','P-R-Met_3']:
    if sample.find('T-R')!=-1:
        reference = 'T-S'
    else:
        reference = 'P-S'
        
    fout.write(
        sample + ',' + reference + ',' +
        bam_dir + sample + '/' + sample + '.gatk.bam,' +
        bam_dir + sample + '/' + sample + '.gatk.bam.bai,' +
        bam_dir + reference + '/' + reference + '.gatk.bam,' +
        bam_dir + reference + '/' + reference + '.gatk.bam.bai\n'
    )

fout.close()
