import pickle
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys
import os

pairs_pth = sys.argv[1]
individual_pth = sys.argv[2]
cell_types = pd.read_csv(sys.argv[3])['cell_types'].to_numpy()
pths = [pairs_pth, individual_pth]
images = []

for k in range(2):
   file_pth = pths[k]
   out = pickle.load(open(file_pth, 'rb'))
   prec_mats = [d['prec_mat'] for d in out]
   c_types = prec_mats[0].shape[0]
   runs_ = list()
   for c in range(c_types):
       tmp = np.array([])
       for run in range(len(prec_mats)):
         tmp = np.append(tmp,prec_mats[run][c][c])
       runs_.append(tmp)
   images.append(runs_)
           
map_idx = {0:'Pairwise TFBMs',1:'Individual TFBMs'}
img_flatten = {}         
for i in range(2):
	tmp = []
	for j in range(len(images[i])):
		tmp += list(images[i][j])
	img_flatten[map_idx[i]] = tmp

category = []

for i in range(len(images[0])):
  category += [cell_types[i]]*len(images[0][0])
           
data = {'Cell types': category,'Pairwise TFBMs': img_flatten['Pairwise TFBMs'],'Individual TFBMs': img_flatten['Individual TFBMs']}

df =  pd.DataFrame(data)

# Melt the dataframe to long format for seaborn
df_melted = pd.melt(df, id_vars=['Cell types'], value_vars=['Pairwise TFBMs', 'Individual TFBMs'], 
                    var_name='Method', value_name='Precision')

# Create the boxplot
plt.figure(figsize=(3, 3))
sns.boxplot(x='Cell types', y='Precision', hue='Method', data=df_melted)

plt.xlabel('Cell types',fontsize=10)
plt.ylabel('Precision',fontsize=10)
dir_path = "../smlr_results"
os.makedirs(dir_path, exist_ok=True)  # Avoids error if directory exists
plt.title('Error Boxplots', fontsize=10)
plt.tight_layout()
plt.savefig('../smlr_results/boxplot_errorbars.pdf',format='pdf',transparent=True,bbox_inches='tight')
