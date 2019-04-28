# pathway interactions
import pandas as pd
import numpy as np
import os

df = pd.read_csv('GSE18123GPL570.csv')
feats = list(df.columns)
feats.remove('condition')
feat_idx = {}
for i in range(len(feats)):
    feat_idx[feats[i]] = i
dirin = 'mapping/'
methods = [file.replace('.csv','') for file in os.listdir(dirin)  if file.endswith(".csv")]
for method in methods:
    print(method)
    connectivity = np.zeros((len(feats),len(feats)))
    for i in range(len(feats)):
        connectivity[i][i] = 1
    df = pd.read_csv(dirin + method + '.csv')
    cols = list(df.columns)
    gene_pathway = dict(zip(df[cols[0]],df[cols[1]]))
    pathway_gene = {}
    for gene in gene_pathway:
        pathways = gene_pathway[gene].split('|')
        for pathway in pathways:
            pathway = pathway.strip()
            if pathway:
                try:
                    pathway_gene[pathway] += [ gene ]
                except:
                    pathway_gene[pathway] = [ gene ]
    for e in pathway_gene:
        genes = list(set(pathway_gene[e]))
        if len(genes)>1:
            for i in range(len(genes)-1):
                i_idx = -1
                try:
                    i_idx = feat_idx[genes[i]]
                except:
                    pass
                if i_idx>=0:
                    for j in range(i+1, len(genes)):
                        j_idx = -1
                        try:
                            j_idx = feat_idx[genes[j]]
                        except:
                            pass
                        if j_idx>=0:
                            connectivity[i_idx][j_idx] = 1
                            connectivity[j_idx][i_idx] = 1
    bf = pd.DataFrame(connectivity)
    bf.to_csv('gene_interaction/' + method.replace('hgncTo','') + '.csv',index=None,header=None)