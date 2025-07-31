import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix
import seaborn as sns
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
from collections import defaultdict
import sys


def shuffle_and_get_indices(x, y):
    # Create a dictionary to store indices for each element in x
    element_indices = defaultdict(list)
    for index, value in enumerate(x):
        element_indices[value].append(index)
    
    # Create the shuffled_x and indices lists
    shuffled_x = []
    indices = []
    for value in y:
        if value in element_indices:
            shuffled_x.extend([value] * len(element_indices[value]))
            indices.extend(element_indices[value])
    
    return shuffled_x, indices

arc_details_pth = sys.argv[1]
model_details_pth = sys.argv[2]
true_labels_pth = sys.argv[3]
motif_pth = sys.argv[4]
c_types = pd.read_csv(sys.argv[5])['cell_types'].to_numpy()

clusters = []
with open(arc_details_pth,'r') as f:
    lines = f.readlines()
    for line in lines:
        clusters.append(int(line.split('\t')[0]))
predicted_classes = np.unique(clusters)
df = pd.read_csv(true_labels_pth)
labels = list(df.iloc[:,0])
true_classes = np.unique(labels)

cm = confusion_matrix(labels, clusters)
cm = cm[~np.all(cm == 0, axis=1)]
cm = cm[:, ~np.all(cm == 0, axis=0)]

# Normalize the confusion matrix by column
column_sums = cm.sum(axis=0, keepdims=True)
normalized_cm = cm / column_sums  # Normalized for coloring
normalized_cm[np.isnan(normalized_cm)] = 0  # Handle division by zero

# Group and reorder predicted clusters
max_true_clusters = np.argmax(cm, axis=0)  # Dominant true cluster for each predicted cluster
percentages = cm / column_sums  # Compute percentage for sorting within groups

# Group indices by dominant true cluster
grouped_indices = {}
for col_idx, true_cluster in enumerate(max_true_clusters):
    if true_cluster not in grouped_indices:
        grouped_indices[true_cluster] = []
    grouped_indices[true_cluster].append(col_idx)

# Within each group, reorder by descending percentage
sorted_indices = []
for true_cluster in sorted(grouped_indices.keys()):  # Sort groups by true cluster order
    cols_sorted = sorted(
        grouped_indices[true_cluster],
        key=lambda col_idx: percentages[true_cluster, col_idx],
        reverse=False,
    )
    sorted_indices.extend(cols_sorted)  # Append columns in correct order

# Reorder the confusion matrix and labels
reordered_cm = cm[:, sorted_indices]
normalized_cm_reordered = normalized_cm[:, sorted_indices]
predicted_labels_sorted = np.array([f'Cluster {i}' for i in np.arange(1,cm.shape[1]+1)[sorted_indices]])

# Custom y-tick labels
true_labels_custom = c_types

# Plot the heatmap
plt.figure(figsize=(8, 6))
cax = plt.imshow(normalized_cm_reordered, interpolation='nearest', cmap='Blues')

# Add a resized colorbar
plt.colorbar(cax, fraction=0.02, pad=0.04)

# Set x and y ticks
plt.xticks(ticks=np.arange(len(predicted_labels_sorted)), labels=predicted_labels_sorted)
plt.yticks(ticks=np.arange(cm.shape[0]), labels=true_labels_custom)

# Add annotations (absolute values)
for i in range(reordered_cm.shape[0]):
    for j in range(reordered_cm.shape[1]):
        plt.text(j, i, f"{reordered_cm[i, j]}", ha="center", va="center", color="black")

plt.xlabel("Predicted Labels",labelpad=15,fontsize=12)
plt.ylabel("True Labels",labelpad=10,fontsize=12)
plt.title("Confusion Matrix")
dir_path = "../nplb_results"
os.makedirs(dir_path, exist_ok=True)
plt.savefig('../nplb_results/confusion_matrix.pdf',format='pdf')

motifs = pd.read_csv(motif_pth)['motifs'].to_numpy()
c_features = {}

imp_features = set()
with open(model_details_pth,'r') as f:
    lines = f.readlines()
    for j, line in enumerate(lines):
        if len(line.split()) > 0 and line.split()[0] == 'Architecture':
            key = int(line.split()[1][0])     
            if key != 1:
                zero_idx = np.argsort(np.array(f_zero))
                one_idx = np.argsort(np.array(f_one))
                c_features[key-1] = [] 
                for i in one_idx[-3:]:
                    imp_features.add(f[i])
                    c_features[key-1].append(motifs[f[i]-1])
                for i in zero_idx[-3:]:
                    imp_features.add(f[i])
                    c_features[key-1].append(motifs[f[i]-1])

            f_zero = []
            f_one = []
            f = []

        if len(line.split()) > 0 and line.split()[0].isdigit():
            f_zero.append(float(line.split()[1][1:]))
            f_one.append(float(line.split()[3][1:]))
            f.append(int(line.split()[0]))
        

        if j == len(lines)-2:
            zero_idx = np.argsort(np.array(f_zero))
            one_idx = np.argsort(np.array(f_one))
            c_features[key] = []
            for i in one_idx[-3:]:
                imp_features.add(f[i])
                c_features[key].append(motifs[f[i]-1])
            for i in zero_idx[-3:]:
                imp_features.add(f[i])
                c_features[key].append(motifs[f[i]-1])

imp_features = list(imp_features)

features_wu_pairs = []
for i in imp_features:
    features_wu_pairs.append(motifs[i-1])
    
str_int = {'C':1, 'A': 0}
with open(arc_details_pth,'r') as f:
    lines = f.readlines()
    data = np.zeros((len(lines),len(imp_features)))
    cluster_labels = np.array([])
    for j, line in enumerate(lines):
        s = line.split('\t')[2]
        c = int(line.split('\t')[0])
        cluster_labels = np.append(cluster_labels,c)
        subset = np.array([str_int[s[i-1]] for i in imp_features])
        data[j] = subset

sorted_labels, indices = shuffle_and_get_indices(cluster_labels, np.arange(1,cm.shape[1]+1)[sorted_indices])
sorted_data = data[indices, :]
print(sorted_labels[:10])

unique_clusters, cluster_starts = np.unique(sorted_labels, return_index=True)
cluster_starts = cluster_starts[sorted_indices]
cluster_ends = np.append(cluster_starts[1:], len(sorted_labels))  # End index for each cluster
cluster_midpoints = (cluster_starts + cluster_ends) // 2

# Create cluster names (e.g., 'Cluster 1', 'Cluster 2', ...)
cluster_names = predicted_labels_sorted
print(cluster_names)

# Find cluster boundaries
boundaries = np.where(np.diff(sorted_labels))[0] + 1

# Plot heatmap
plt.figure(figsize=(8, 6))

cmap = sns.color_palette("viridis", as_cmap=True)

# Create a discrete normalization for binary data
bounds = [0, 0.5, 1]  # Boundaries for 0 and 1
norm = mcolors.BoundaryNorm(bounds, ncolors=256)

sns.heatmap(sorted_data, cmap=cmap, norm=norm, cbar_kws={'ticks': [0, 1],"shrink": 0.5}, xticklabels=features_wu_pairs, yticklabels=False)

plt.yticks(cluster_midpoints + 0.5, cluster_names, rotation=0,fontsize=7)

plt.xticks(rotation=90,fontsize=6)

for b in boundaries:
    plt.axhline(b, color='white', linestyle='--', linewidth=1.5)

plt.title("Heatmap with Cluster Boundaries")
plt.xlabel("Motifs")
plt.ylabel("Cells")

plt.savefig('../nplb_results/heatmap.pdf',format='pdf')
