from sklearn.model_selection import GridSearchCV
import numpy as np
import pandas as pd
import sys
import pickle
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import confusion_matrix
import os
import re

data_pth = sys.argv[1]
labels_pth = sys.argv[2]

data = pd.read_csv(data_pth)
X = data.to_numpy()

labels = pd.read_csv(labels_pth)
Y = labels['labels'].to_numpy()

clf = LogisticRegression(
            solver="saga",
            penalty=None,
            random_state=42,
        )

clf.fit(X,Y)
Y_predict = clf.predict(X)

conf_mat = confusion_matrix(Y, Y_predict)

out = dict(conf_mat=conf_mat,clf=clf,motifs=np.array(data.columns))
dir_path = "../smlr_results"
os.makedirs(dir_path, exist_ok=True)  
match = re.search(r'([^/]+)\.csv$', data_pth)
if match:
    fname = match.group(1)
print(fname)
out_pkl_file = "../smlr_results/{fname}_feature_weights_full.pkl".format(fname=fname)
print(out_pkl_file)

with open(out_pkl_file, 'wb') as file:  
    pickle.dump(out, file)

