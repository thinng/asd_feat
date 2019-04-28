import pandas as pd
import numpy as np
ls = ['co_expression','DO', 'GObp','GOcc', 'GOmf', 'HPO','PPI','pubmed']
total = None
for i in range(len(ls)):
    print(i,ls[i])
    df = pd.read_csv('gene_interaction/' + ls[i]+'.csv', header=None)
    if i==0:
        total = df.values
    else:
        total = np.sum([total, df.values], axis=0)
x = total/len(ls)
x = (x>0.5).astype(int)
bf = pd.DataFrame(x)
bf.to_csv('gene_interaction/joint_interaction.csv',index=None,header=None)