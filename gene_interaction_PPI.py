# PPI interaction
import pandas as pd
import numpy as np

df = pd.read_csv('GSE18123GPL570.csv')
feats = list(df.columns)
feats.remove('condition')

df = pd.read_csv('mapping/hippie_current.txt', header=None, sep='\t')
df.columns = ['p1','v1','p2','v2','score','desc']
df['p1'] = df['p1'].apply(lambda x: str(x).replace('_HUMAN',''))
df['p2'] = df['p2'].apply(lambda x: str(x).replace('_HUMAN',''))
df = df[ ['p1','p2','score'] ]
connectivity = np.zeros((len(feats),len(feats)))
for i in range(len(feats)):
    connectivity[i][i] = 1
for index, row in df.iterrows():
    print , row["p2"]
    i_idx = -1
    try:
        i_idx = feat_idx[row["p1"]]
    except:
        pass
    if i_idx>=0:
        j_idx = -1
        try:
            j_idx = feat_idx[row["p2"]]
        except:
            pass
        if j_idx>=0:
            connectivity[i_idx][j_idx] = 1
            connectivity[j_idx][i_idx] = 1
bf = pd.DataFrame(connectivity)
bf.to_csv('gene_interaction/PPI.csv',index=None,header=None)