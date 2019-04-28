# PubMed interaction
#import Python libaries needed 
from gensim.models import Word2Vec
import pandas as pd
import numpy as np

def get_neighbor(gene):
    ret = ''
    try:
        x = model.most_similar(gene, topn=10)
    except:
        pass
    else:
        similar_list = [ e[0] for e in x]
        tmp = list(set(gene_list).intersection(set(similar_list)))
        tmp = [genename_emb_dataset[e] for e in tmp if e in genename_emb_dataset]
        if len(tmp)>0:
            ret = ','.join(map(str,tmp))
    return ret


df = pd.read_csv('GSE18123GPL570.csv')
feats = list(df.columns)
feats.remove('condition')

# As PubMed data is in lowercase for learning embeddings, 
# we map lower case gene names to the name in the dataset 
genename_emb_dataset = {}
for line in open('pubmed/gene_name_ref.txt'):
    cols = line.strip().split()
    if len(cols)==2:
        if cols[0] in feats:
            genename_emb_dataset[ cols[1] ] = cols[0]
            
model = Word2Vec.load('pubmed/gene_embedding.embed')
        
gene_list = list(genename_emb_dataset.keys()) 
df = pd.DataFrame(gene_list)
df.columns = ['embed_name']
df['dataset_name'] = df['embed_name'].apply(lambda x: genename_emb_dataset[x])
df['neighbors'] = df['embed_name'].apply(get_neighbor)
df = df[ df['neighbors'] !='' ]

feat_idx = {val:idx for idx, val in enumerate(feats)}
    
connectivity = np.zeros((len(feats),len(feats)))
for i in range(len(feats)):
    connectivity[i][i] = 1
    
for index, row in df.iterrows():
    i_idx = feat_idx[row["dataset_name"]]
    neighbor = row['neighbors']
    genes = neighbor.split(',')
    for gene in genes:
        j_idx = feat_idx[gene]
        connectivity[i_idx][j_idx] = 1
        connectivity[j_idx][i_idx] = 1
bf = pd.DataFrame(connectivity)
bf.to_csv('gene_interaction/pubmed.csv',index=None,header=None)
