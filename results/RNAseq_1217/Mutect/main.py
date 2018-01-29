import vcf
import pandas as pd
import seaborn as sns

fin = vcf.Reader(open('mutect.vcf','r'))
data = []
for v in fin:
    tmp = []
    if len(v.ALT)>1:
        continue

    n = 0
    nt = 0
    nn = 0
    for s in v.samples:
        if s.called:
            n+=1
            if s.sample in ['R1-HCY1','R1-HCY2','R1-HCY3','R1-MET1','R1-MET2','R1-MET3','R2-HCY1','R2-HCY2','R2-HCY3','R2-MET1','R2-MET2','R2-MET3']:
                nt+=1
            else:
                nn+=1
            try:
                tmp.append(float(s.data.AF))
            except:
                import pdb; pdb.set_trace()
        else:
            tmp.append(0.0)
    if n>4 and n < 12  and nt > 2 and nn > 2:
        data.append(tmp)
data = pd.DataFrame(data)
data.columns = [s.sample for s in v.samples]
data.to_csv('test.csv',index=None)

sns.set(font_scale=0.5)
g = sns.clustermap(data)
g.savefig('test.png')
