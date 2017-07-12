import pandas
import pyensembl

def parse_result_file(filename):
    mapping = []
    raw = []
    bed = []
    index=0

    data = open(filename,'r').read().split('\n')
    mapping.append(data[index].replace('<mapping>','').split(','))
    index+=1
    while data[index].find('</mapping><raw>')==-1:
        mapping.append(data[index].split(','))
        index+=1
    raw.append(data[index].replace('</mapping><raw>','').split('\t'))
    index+=1
    while data[index].find('</raw><bed>')==-1:
        raw.append(data[index].split('\t'))
        index+=1
    index+=1

    while data[index].find('</bed>')==-1:
        bed.append(data[index].split('\t'))
        index+=1

    mapping = pandas.DataFrame(mapping,columns=['gene_id1','gene_id2','chrom','start','end','gene_id3'])
    raw = pandas.DataFrame(raw,columns=['chrom','start','end','pvalue','padj','num_differentially_expressed_genes','total_genes_in_region'])
    bed = pandas.DataFrame(bed,columns=['chrom','start','end','pval'])

    return mapping,raw,bed


def prep_datat():
    deseq = pandas.read_csv("/Users/charlesmurphy/Desktop/Research/0617_yuxiang/results/RNAseq/DESeq2/parental_vs_resistant_ensembl/parental_vs_resistant_results.csv")
    deseq2 = deseq[(deseq['padj']<=0.01) & (abs(deseq['log2FoldChange']) >= 1)]
    deseq2['Unnamed: 0'].to_csv('parental_vs_resistant.txt',index=False,sep='\t')

    deseq2 = deseq[(deseq['padj']<=0.01) & (deseq['log2FoldChange'] >= 1)]
    deseq2['Unnamed: 0'].to_csv('parental_vs_resistant_up.txt',index=False,sep='\t')

    deseq2 = deseq[(deseq['padj']<=0.01) & (deseq['log2FoldChange']<= -1)]
    deseq2['Unnamed: 0'].to_csv('parental_vs_resistant_down.txt',index=False,sep='\t')


    deseq = pandas.read_csv("/Users/charlesmurphy/Desktop/Research/0617_yuxiang/results/RNAseq/DESeq2/parental_vs_resistant_Met-Hcy+_ensembl/parental_vs_resistant_Met-Hcy+_results.csv")
    deseq = deseq[(deseq['padj']<=0.01) & (abs(deseq['log2FoldChange']) >= 1)]
    deseq2['Unnamed: 0'].to_csv('parental_vs_resistant_Met-Hcy+.txt',index=False,sep='\t')

    deseq2 = deseq[(deseq['padj']<=0.01) & (deseq['log2FoldChange'] >= 1)]
    deseq2['Unnamed: 0'].to_csv('parental_vs_resistant_Met-Hcy+_up.txt',index=False,sep='\t')

    deseq2 = deseq[(deseq['padj']<=0.01) & (deseq['log2FoldChange']<= -1)]
    deseq2['Unnamed: 0'].to_csv('parental_vs_resistant_Met-Hcy+_down.txt',index=False,sep='\t')


    deseq = pandas.read_csv("/Users/charlesmurphy/Desktop/Research/0617_yuxiang/results/RNAseq/DESeq2/parental_vs_resistant_Met+Hcy-_ensembl/parental_vs_resistant_Met+Hcy-_results.csv")
    deseq = deseq[(deseq['padj']<=0.01) & (abs(deseq['log2FoldChange']) >= 1)]
    deseq2['Unnamed: 0'].to_csv('parental_vs_resistant_Met+Hcy-.txt',index=False,sep='\t')

    deseq2 = deseq[(deseq['padj']<=0.01) & (deseq['log2FoldChange'] >= 1)]
    deseq2['Unnamed: 0'].to_csv('parental_vs_resistant_Met+Hcy-_up.txt',index=False,sep='\t')

    deseq2 = deseq[(deseq['padj']<=0.01) & (deseq['log2FoldChange']<= -1)]
    deseq2['Unnamed: 0'].to_csv('parental_vs_resistant_Met+Hcy-_down.txt',index=False,sep='\t')

def getOverlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))

def get_overlapping_genes(raw,mapping,genes):
    overlapping_genes = []
    mapping['chrom']=mapping['chrom'].str.upper()
    for i in raw.index:
        region_start = int(raw.loc[i,'start'])
        region_end = int(raw.loc[i,'end'])
        chrom = raw.loc[i,'chrom']
        mapping_sub = mapping[mapping['chrom']==chrom]

        subset = []
        for j in mapping_sub.index:
            j_start = int(mapping_sub.loc[j,'start'])
            j_end = int(mapping_sub.loc[j,'end'])
            if getOverlap([region_start,region_end],[j_start,j_end])>0:
                subset.append(True)
            else:
                subset.append(False)

        mapping_sub = mapping_sub[subset]['gene_id1'].tolist()
        mapping_sub = map(lambda x: x.upper(),mapping_sub)
        mapping_sub = ','.join(genes.loc[mapping_sub,'symbols'].tolist())
        overlapping_genes.append(mapping_sub)

    raw['overlapping_differentially_expressed_genes'] = overlapping_genes
    raw['padj'] = raw['padj'].astype(float)
    raw = raw.sort_values(['padj'])
    return raw

def process_results():

    db = pyensembl.EnsemblRelease('75','human')

    for directory in ['parental_vs_resistant','parental_vs_resistant_Met-Hcy+','parental_vs_resistant_Met+Hcy-']:

        genes = pandas.read_table(directory + '.txt',header=None)
        symbols = []
        for ii in genes[0].tolist():
            try:
                symbols.append(db.gene_by_id(ii).gene_name)
            except ValueError:
                symbols.append(ii)
        genes['symbols'] = symbols
        genes.index = genes[0]

        mapping,raw,bed = parse_result_file(
            './' + directory + '/' + directory + '_results.txt'
        )
        raw = get_overlapping_genes(raw,mapping,genes)
        raw.to_excel('./' + directory + '/' + directory + '_results.xlsx',index=False)


        genes = pandas.read_table(directory + '_up.txt',header=None)
        symbols = []
        for ii in genes[0].tolist():
            try:
                symbols.append(db.gene_by_id(ii).gene_name)
            except ValueError:
                symbols.append(ii)
        genes['symbols'] = symbols
        genes.index = genes[0]

        mapping,raw,bed = parse_result_file(
            './' + directory + '/' + directory + '_results_up.txt'
        )
        raw = get_overlapping_genes(raw,mapping,genes)
        raw.to_excel('./' + directory + '/' + directory + '_results_up.xlsx',index=False)


        genes = pandas.read_table(directory + '_down.txt',header=None)
        symbols = []
        for ii in genes[0].tolist():
            try:
                symbols.append(db.gene_by_id(ii).gene_name)
            except ValueError:
                symbols.append(ii)
        genes['symbols'] = symbols
        genes.index = genes[0]

        mapping,raw,bed = parse_result_file(
            './' + directory + '/' + directory + '_results_down.txt'
        )
        raw = get_overlapping_genes(raw,mapping,genes)
        raw.to_excel('./' + directory + '/' + directory + '_results_down.xlsx',index=False)

#prep_datat()
process_results()
