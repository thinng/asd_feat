import pandas as pd
import numpy as np
from skfeature.function.information_theoretical_based import CMIM,ICAP,JMI
from skfeature.function.similarity_based import fisher_score,reliefF,trace_ratio
from skfeature.function.statistical_based import chi_square, f_score, gini_index
from skfeature.function.sparse_learning_based import ll_l21,ls_l21,RFS
from skfeature.utility.sparse_learning import *
from sklearn import metrics
import sys

def normalized(x):
    return (x-min(x))/(max(x)-min(x))

df = pd.read_csv('GSE18123GPL570.csv')
y = list(map(int,list(df['condition']))) 
y_train = np.asarray(y)

feats = list(df.columns) 
feats.remove('condition')
df = df[feats]
max_feat = len(feats)
X_train = df.as_matrix()

sel = ['chi_square','cmim', 'f_score','fisher_score','gini_index','icap',
       'jmi','ll_l21','ls_l21','reliefF','rfs','trace_ratio', 'SFARI','AUC'] [int(sys.argv[1])]
print('running feature selection method:', sel)
if sel=='chi_square':
    score = chi_square.chi_square(X_train, y_train)
    idx = chi_square.feature_ranking(score)
elif sel=='cmim':
    idx,score, MIfy = CMIM.cmim(X_train, y_train)
elif sel=='f_score':
    score = f_score.f_score(X_train, y_train)
    idx = f_score.feature_ranking(score)
elif sel=='fisher_score':
    score = fisher_score.fisher_score(X_train, y_train)
    idx = fisher_score.feature_ranking(score)
elif sel=='gini_index':
    score = gini_index.gini_index(X_train, y_train)
    idx = gini_index.feature_ranking(score)
elif sel=='icap':
    idx,score, MIfy = ICAP.icap(X_train, y_train)
elif sel=='jmi':
    idx,score, MIfy = JMI.jmi(X_train, y_train)
elif sel=='ll_l21':
    score, obj, value_gamma = ll_l21.proximal_gradient_descent(X_train, construct_label_matrix_pan(y_train), 0.1, verbose=False)
    idx = feature_ranking(score)
elif sel=='ls_l21':
    score, obj, value_gamma = ls_l21.proximal_gradient_descent(X_train, construct_label_matrix_pan(y_train), 0.1, verbose=False)
    idx = feature_ranking(score)
elif sel=='reliefF':
    score = reliefF.reliefF(X_train, y_train)
    idx = reliefF.feature_ranking(score)
elif sel=='rfs':
    score = RFS.rfs(X_train, construct_label_matrix(y_train), gamma=0.1)
    idx = feature_ranking(score)
elif sel=='trace_ratio':
    idx, score, subset_score = trace_ratio.trace_ratio(X_train, y_train,max_feat, style='fisher') 
elif sel=='SFARI':
    bf = pd.read_csv('mapping/SFARI-Gene_genes_export09-06-2018.csv')
    gene_reports = dict(zip(bf['gene-symbol'], bf['number-of-reports']))
    sorted_ = sorted(gene_reports.items(), key=lambda kv: kv[1], reverse=True)
    gene_sorted = [e[0] for e in sorted_]    
    sfari_idx = [feats.index(e) for e in gene_sorted if e in feats]
    non_sfari_idx = [i for i in range(len(feats)) if not i in sfari_idx]
    idx = sfari_idx + non_sfari_idx    
elif sel=='AUC':
    idx_auc = {}
    for i in range(len(feats)):
        y_predict = normalized(X_train[:,i])
        fp, tp, _ = metrics.roc_curve(y_train, y_predict)
        auc = metrics.auc(fp, tp)
        idx_auc[i] = auc
    sorted_ = sorted(idx_auc.items(), key=lambda kv: kv[1], reverse=True)
    idx = [e[0] for e in sorted_]    
else:
    print('wrong feature selection!!!')
with open('feat_ranking/' + sel + '.txt', 'w') as f:
    f.write('\n'.join(map(str,idx)))
print('Done feature selection method:', sel)