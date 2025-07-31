import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
import re
import sys
import os

pairs_pth = sys.argv[1]
individual_pth = sys.argv[2]
pths = [pairs_pth, individual_pth]

out = pickle.load(open(pths[0], 'rb'))
c = out['clf'].coef_ 
n_c_types = c.shape[0]

fig, axes = plt.subplots(n_c_types, 2, figsize=(8, 10))  # Create a 2x3 grid of subplots

for i in range(2):  # Outer loop (2 iterations)
   file_pth = pths[i]
   out = pickle.load(open(file_pth, 'rb'))
   motifs = out['motifs']
   c = out['clf'].coef_ 
   for j in range(3):  # Inner loop (3 iterations)
      ax = axes[j, i]
      t = c[j,]
      indices = np.argsort(t)
      x = motifs[indices]
      y = t[indices]
  
      x_indices = range(len(x))

      label_indices = list(range(5)) + list(range(len(x) - 5, len(x)))

      ax.scatter(x_indices,y,color='blue', alpha=0.7)

      n = len(indices)

      min_y_gap = 0.10  

      adjusted_positions = []

      mots = []

      for k in label_indices:
    # Check and adjust for vertical overlap
         adjusted_y = y[k]
         if len(adjusted_positions)>0:
            if (adjusted_y <= adjusted_positions[-1]) or ((adjusted_y - adjusted_positions[-1]) < min_y_gap):
               adjusted_y = adjusted_positions[-1] + min_y_gap  # Shift the annotation upwards
         print(adjusted_y)
         adjusted_positions.append(adjusted_y)  # Store adjusted position

         ax.annotate(
            "{x_val}, {y_val:.2f}".format(x_val=x[k],y_val=y[k]),       # Annotation text
            xy=(k, adjusted_y),                # Coordinates of the point
            textcoords="offset points",  # Position relative to the point
            xytext=(0,0),  # Dynamically offset text
            fontsize=7,
            fontweight='semibold',
            arrowprops=dict(arrowstyle="->", color='gray', lw=0.5)
         )


      ax.set_xlabel('Motifs',fontsize=13)
      ax.set_ylabel('Feature importance',fontsize=13)
      ax.set_title('MLR Feature Weights',fontsize=15)

      ax.spines['top'].set_visible(False)
      ax.spines['right'].set_visible(False)
   
plt.tight_layout()
plt.subplots_adjust(hspace=0.5, wspace=0.3)
dir_path = "../smlr_results"
os.makedirs(dir_path, exist_ok=True)  # Avoids error if directory exists
plt.savefig('../smlr_results/feature_weights.pdf',format = 'pdf',bbox_inches='tight')

# plt.savefig('feature_importance_smlr_furlong_same_noPreproc_C_0.1_10_12hr.png',transparent=True)

