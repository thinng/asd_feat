# co_expression interactions
import pandas as pd
import numpy as np

df = pd.read_csv('GSE18123GPL570.csv')
feats = list(df.columns) 
feats.remove('condition')
df = df[feats]
df.to_csv('gene_interaction/expression.csv',index=None,header=None)
vv = df.corr().values
connectivity = np.zeros(vv.shape)
for i in range(len(feats)-1):
    for j in range(i+1,len(feats)):
        w = abs(vv[i][j])
        if w>=0.5:
            connectivity[i][j] = 1
            connectivity[j][i] = 1
for i in range(len(feats)):
    connectivity[i][i] = 1
bf = pd.DataFrame(connectivity)
bf.to_csv('gene_interaction/co_expression.csv',index=None,header=None)