from sklearn.model_selection import GridSearchCV
import re
import numpy as np
import pandas as pd
import sys
import pickle
from sklearn.metrics import classification_report, confusion_matrix, f1_score
from sklearn.model_selection import StratifiedKFold
from joblib import Parallel, delayed
from sklearn.linear_model import LogisticRegression
import os

data_pth = sys.argv[1]
labels_pth = sys.argv[2]

data = pd.read_csv(data_pth)
labels = pd.read_csv(labels_pth)

X = data.to_numpy()
Y = labels['labels'].to_numpy()

clf = LogisticRegression(
            solver="saga",
            penalty=None,
            random_state=42,
        )
        
n_labels = len(np.unique(Y))

def eachRun(i):
    skf = StratifiedKFold(n_splits=10,shuffle=True,random_state=i)
    res = np.zeros((n_labels,n_labels))
    for train_index, test_index in skf.split(X, Y):
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = Y[train_index], Y[test_index]

        clf.fit(X_train, y_train)
        Y_predict = clf.predict(X_test)

        #1 = f1_score(y_test, Y_predict, average='weighted')
        conf_mat = confusion_matrix(y_test, Y_predict)
        res = np.add(res,conf_mat)
    prec = (res.T/res.sum(axis=1)).T
        
    return dict(conf_mat=res,prec_mat=prec)
    
n_runs = 10

out = Parallel(n_jobs=n_runs, verbose=100)(
      delayed(eachRun)(i) for i in range(1,n_runs+1))

dir_path = "../smlr_results"
os.makedirs(dir_path, exist_ok=True)
match = re.search(r'([^/]+)\.csv$', data_pth)
if match:
    fname = match.group(1)
out_pkl_file = "../smlr_results/{fname}_precision_full_kfold_multiple_runs.pkl".format(fname=fname)

with open(out_pkl_file, 'wb') as file:  
    pickle.dump(out, file)

