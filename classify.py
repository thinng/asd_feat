import pandas as pd
import numpy as np
from sklearn import metrics
from sklearn.cross_decomposition import PLSRegression
from sklearn import model_selection
from sklearn.preprocessing import scale 
import os
            
def evaluate(y_true, y_pred):
    fp, tp, _ = metrics.roc_curve(y_true, y_pred)
    auc = metrics.auc(fp, tp)
    ret = [metrics.accuracy_score(y_true, y_pred)] + [auc] + [metrics.f1_score(y_true, y_pred)] + [metrics.matthews_corrcoef(y_true, y_pred)]
    return ret

datasets = ['GSE18123GPL570','GSE18123GPL6244']
df = pd.read_csv(datasets[0] + '.csv')
y = list(map(int,list(df['condition'])))
y = [-1 if e==0 else e for e in y]
y_train = np.asarray(y)
feats = list(df.columns) 
feats.remove('condition')
max_feat = len(feats)
df = df[feats]
X_train = df.as_matrix()

df = pd.read_csv(datasets[1] + '.csv')
y = list(map(int,list(df['condition']))) 
y = [-1 if e==0 else e for e in y]
y_test = np.asarray(y)
df = df[feats]
X_test = df.as_matrix()
X_train, X_test = scale(X_train), scale(X_test) 

cf = pd.read_csv('mapping/SFARI-Gene_genes_export09-06-2018.csv')
sfari = list(cf['gene-symbol'])
sfari_idx = [feats.index(e) for e in sfari if e in feats]


kf_10 = model_selection.KFold(n_splits=10, shuffle=True, random_state=1)

dirin = 'feat_ranking/'
methods = []
for file in os.listdir(dirin):
    if file.endswith(".txt"):
        methods += [ file.replace('.txt','') ]

print('number of methods:', len(methods))

result = []
for i in range(len(methods)):
    method = methods[i]
    print('runnning at ', i, method)
    idx = [int(e) for e in open(dirin + method + '.txt') if e.strip()]
    for jj in range(50):
        num_feat =  (jj+1)*10 
        train_X = X_train[:,idx[0:num_feat]]
        test_X = X_test[:,idx[0:num_feat]]
#         learning out the number of components providing best performance from training data
        mse = []
        for k in range(10):
            pls = PLSRegression(n_components=k+1)
            score = model_selection.cross_val_score(pls, train_X, y_train, cv=kf_10, scoring='neg_mean_squared_error').mean()
            mse.append(score)
        n_component = np.argmax(mse) + 1
        clf = PLSRegression(n_components=n_component)
        clf.fit(train_X, y_train)
        y_predict = clf.predict( test_X )
        y_predict = [1 if e>=0 else -1 for e in y_predict]
        ret = evaluate(y_test, y_predict)
        ret = [round(e,3) for e in ret]
        in_safari = len(set(sfari_idx).intersection(set(idx[0:num_feat])))
        result += [ [method, num_feat, in_safari] + ret  ]
    # if i==1:
        # break

df = pd.DataFrame(result)
df.columns = ['method','number_feature','in_safari', 'accuracy','auc','f1','mcc']
df.to_csv('result.csv', index=None)