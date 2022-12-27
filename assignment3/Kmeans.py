import pandas as pd
import numpy as np
from sklearn import preprocessing
from sklearn.cluster import KMeans
from math import sqrt
import scipy

file = pd.read_csv("/Users/sarapour/Desktop/Biase_2014.csv", sep=",", index_col=0).T

scaler = preprocessing.StandardScaler()
scaled_data = scaler.fit_transform(file)
file_scaled = pd.DataFrame(scaled_data)
file_scaled.columns = file.columns
file_scaled['sample_id'] = file.index

def getRandomCentroids(array, num_clus, random_state):
    x = array
    x = x[x.columns.difference(['sample_id'])]
    x = np.array(x)
    np.random.seed(random_state)
    random_idx = np.random.permutation(x.shape[0])
    centroids = x[random_idx[:num_clus]]
    return centroids

# Function: Should Stop
# -------------
# Returns True or False if k-means is done. K-means terminates either
# because it has run a maximum number of iterations OR the centroids
# stop changing.
MAX_ITERATIONS = 20
def shouldStop(oldCentroids, centroids, iterations):  # Whether to stop
    if iterations > MAX_ITERATIONS:
        return True
    return np.array_equal(oldCentroids, centroids)

# Function: Get Centroids
# -------------
# Returns k random centroids, each of dimension n.
def getCentroids(dataSet, k):
    # Each centroid is the geometric mean of the points that
    # have that centroid's label. Important: If a centroid is empty (no points have
    # that centroid's label) you should randomly re-initialize it.
    new_centers = np.array([dataSet[dataSet['labels'] == i].mean(0)
                            for i in range(k)])
    return new_centers

def data_to_centroid_dist(data, centroid):
    # Distance to all points in B, for each point in A.
    dist_mat = scipy.spatial.distance.cdist(data, centroid, 'euclidean')
    # Indexes to minimum distances.
    min_dist_idx = np.argmin(dist_mat, axis=1)
    # Store only the minimum distances for each point in A, to a point in B.
    min_dists = [dist_mat[i][min_d_idx] for i, min_d_idx in enumerate(min_dist_idx)]
    return min_dist_idx, min_dists


# x is the data set, k is the type, the number of cycles up to maxIt
def kmeans(mat, k, rnd):
    # Initialize centroids randomly
    centroids = getRandomCentroids(mat, 4, rnd)
    # Randomly assign labels to initial centorid
    # Initialize book keeping vars.
    iterations = 0
    oldCentroids = None
    # Run the main k-means algorithm
    while not shouldStop(oldCentroids, centroids, iterations):
        # Save old centroids for convergence test. Book keeping.
        oldCentroids = np.copy(centroids)  # Copy an array
        iterations += 1
        # Assign labels to each datapoint based on centroids
        labels, dists = data_to_centroid_dist(mat[mat.columns.difference(['sample_id'])],
                                              centroids)
        mat['labels'] = labels
        # Assign centroids based on datapoint labels
        centroids = getCentroids(mat, 4)[:, :-1]
        mat.drop(['labels'], axis = 1, inplace = True)
    # We can get the labels too by calling getLabels(dataSet, centroids)
    return mat, labels


def mann_whit_u(g):
    p_values = np.zeros(((g.shape[1]-2), 4))
    for i in range(4):
        clust_i = g.loc[g[0] == i]
        others_i = g.loc[g[0] != i]
        for j in range((g.shape[1]-2)):
            gene_j_clust = clust_i.iloc[:, j]
            gene_j_others = others_i.iloc[:, j]
            p_val = scipy.stats.mannwhitneyu(gene_j_clust, gene_j_others)[1]
            p_values[j, i] = p_val
    return p_values

#################################CHECK WITH KMEANS##############################
run = kmeans(file_scaled, 4, 9)
sampl = pd.DataFrame(run[0])
labels = pd.DataFrame(run[1])
g = pd.concat([sampl, labels], axis = 1)
clust = g.loc[g[0] == 3]#.drop(['sample_id', 0], axis = 1).T
others = g.loc[g[0] != 3]#.drop(['sample_id', 0], axis = 1).T
p_values = mann_whit_u(g)
# matrix where rows are genes and columns are clust_id vs non_clust_id (e.g. 0 is clust 0 vs 1,2,3)
p_val_df = pd.DataFrame(p_values)
p_val_df.index = g.columns[:-2]
# bonferroni correct for 4 tests: 0.05/4 = 0.0125
# clust 0<- 2 cells
clust_0 = pd.DataFrame(p_val_df.iloc[:, 0])
clust_0.columns = ['p_value']
clust_0 = clust_0[clust_0['p_value'] < 0.0125]
# clust 1 <- 2 cell
clust_1 = pd.DataFrame(p_val_df.iloc[:, 1])
clust_1.columns = ['p_value']
clust_1 = clust_1[clust_1['p_value'] < 0.0125]
# clust 2 <- a couple 4 cell
clust_2 = pd.DataFrame(p_val_df.iloc[:, 2])
clust_2.columns = ['p_value']
clust_2 = clust_2[clust_2['p_value'] < 0.0125]
# clust 3 <- mostly 4 cell
clust_3 = pd.DataFrame(p_val_df.iloc[:, 3])
clust_3.columns = ['p_value']
clust_3 = clust_3[clust_3['p_value'] < 0.0125]

import matplotlib.pyplot as plt
from matplotlib_venn import venn3

set0 = set(clust_0.index)
set1 = set(clust_1.index)
set2 = set(clust_2.index)
set3 = set(clust_3.index)

# len(set0.intersection(set1)) = 299, small_4cell x large_2 cell
# len(set1.intersection(set2)) = 2628, large_2cell x large_4 cell
# len(set2.intersection(set3)) = 505, large_2cell x small_2 cell
# len(set3.intersection(set0)) = 54, small_2 cell x small_4 cell
# len(set3.intersection(set1)) = 380, large_2 cell x small_2 cell
# len(set2.intersection(set0)) = 239, large_4cell x small_4 cell

# for cluster 3, 6812*4 matrix then do cutoff and for each cluster you get list of genes that is signficant. mse/sse for the best result orsay this is mostly 4cell vs 2cell


# #################################

