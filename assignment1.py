

from utilities import import_16s_sequences, print_distance_matrix, heatmap
from utilities import distance_matrix_to_coordinates_MDS, mds_scatterplot
from utilities import silhouette_plot, silhouette_score
from collections import defaultdict
import numpy as np
import warnings
warnings.filterwarnings('ignore')

#### OPENING

### import genes
ribosome_genes = import_16s_sequences()


### Section talking about classes


### provided code block
for gene in ribosome_genes:
    print('{:<10}{:<30}{:<10}{:<10}'.format(gene.accession, gene.name[:27], gene.length, gene.sequence[:8] + '...'))



##############################################
####### SECTION 1: PAIRWISE DISTANCES ########
##############################################



####### TASK1: kmer distance #######


# The k-mer distance between two sequences is the fraction of k-mers that are unique to either sequence
# Write a function that calculates a dictionary of k-mers (for k = any number) and their counts for each gene sequence.
def create_kmer_dictionary(seq, k):
    kmers = defaultdict(int)

    for i in range(len(seq) - k):
        kmer = seq[i:i + k]
        kmers[kmer] += 1
    
    return kmers

# Write a function that accepts two kmer occurance dictionaries, and returns the kmer distance
def calculate_total_unique_kmers(kmers1, kmers2):
    unique_kmers = 0
    for k in kmers1:
        if k not in kmers2:
            unique_kmers += 1

    for k in kmers2:
        if k not in kmers1:
            unique_kmers += 1

    return unique_kmers


# write a function that uses the above two functions to calcualte the kmer distance of two sequences
def kmer_distance(seq1, seq2, k):
    kmer_dict1 = create_kmer_dictionary(seq1, k)
    kmer_dict2 = create_kmer_dictionary(seq2, k)
    distance = calculate_total_unique_kmers(kmer_dict1, kmer_dict2)
    return distance


# lets check our function. We can use the first two entries in the 'ribosome_genes' list.
# if implemented correctly, the following should return 56
distance = kmer_distance(ribosome_genes[0].sequence, ribosome_genes[1].sequence, 15)
print(distance)





####### TASK2: smith-waterman alignment #######
# Another way to compare the similarity of two sequences is through alignment.
# the alignment score of two sequences will be high when they are similar, and low when they are distinct. 


# write a function that initialises and returns a scoregrid for two sequences
# keep in mind the matrix must be 1 element larger than the sequence lengths, 
# and that the indel scores for the first row and column need to be filled in
def init_scoregrid(seq1, seq2, indel_score=-4):
    xdim = len(seq1) + 1
    ydim = len(seq2) + 1
    scoregrid = np.zeros((ydim, xdim), np.int)
    scoregrid[0, :] = list(range(0, indel_score * xdim, indel_score))
    scoregrid[:, 0] = list(range(0, indel_score * ydim, indel_score))
    return scoregrid


# lets do a sanity check that the grid has been initialised properly. 
# The folling should print the initialised scoregrid
print(init_scoregrid('hello', 'kittycat'))



# write a function that calculates the initialised scoregrid. accepts two sequences and a scoregrid. 
def calculate_scoregrid(seq1, seq2, scoregrid, match_score=1, mismatch_score=-4, indel_score=-4):
    xdim = len(seq1) + 1
    ydim = len(seq2) + 1

    for x in range(1, xdim):
        for y in range(1, ydim):
            diagonal_score = match_score if seq1[x - 1] == seq2[y - 1] else mismatch_score
            score = max(scoregrid[y - 1, x] + indel_score,
                        scoregrid[y, x - 1] + indel_score,
                        scoregrid[y - 1, x - 1] + diagonal_score)
            scoregrid[y, x] = score

    return scoregrid


# lets do another sanity check 
# The folling should print a calculated scoregrid, with the following in the bottom right:
#   0  0 
#  -4  1
scoregrid = init_scoregrid('hello', 'helllo')
print(calculate_scoregrid('hello', 'helllo', scoregrid))


# write a function that returns the alignment score (smith-waterman) from a scoregrid
def report_alignment_score(scoregrid):
    return np.max(scoregrid)


# Final sanity check. Should return 4.
scoregrid = init_scoregrid('hello', 'helllo')
calculated_scoregrid = calculate_scoregrid('hello', 'helllo', scoregrid)
print(report_alignment_score(calculated_scoregrid))


# write a function that uses the above three functions to calculate the alignment score of two sequences
def smith_waterman(seq1, seq2):
    scoregrid = init_scoregrid(seq1, seq2)
    calculated_scoregrid = calculate_scoregrid(seq1, seq2, scoregrid)
    alignment_score = report_alignment_score(calculated_scoregrid)
    return alignment_score


# the following should print 4
print(smith_waterman('hello', 'helllo'))


####### TASK3: pairwise distances #######

# We have now written two functions which can calculate the distance of two sequences. 
# We can calculate the kmer distance, and the smith-waterman alignment score. 
# lets use these two methods to calculate the pairwise distance of our genes. 

# write a function which initialises and returns a distance matrix from a list of genes. 
# use a dictionary, where the keys are the gene.accession, and the values are all 0 for now. 

def init_distance_matrix(genes):
    distance_matrix = defaultdict(dict)
    for gene1 in genes:
        for gene2 in genes:
            distance_matrix[gene1.accession][gene2.accession] = 0
    
    return distance_matrix


# Lets print the distance matrix to make sure it worked. 
#distance_matrix = init_distance_matrix(ribosome_genes)
#print_distance_matrix(distance_matrix)


# Time to fill in the matrix with distances. 
# write a function which calculates the pairwise distance of genes using kmer distance.
# you will need to call the 'kmer_distance' function you have written above. 
def calculate_kmer_distance_matrix(genes, matrix, k):
    for gene1 in genes:
        for gene2 in genes:
            matrix[gene1.accession][gene2.accession] = kmer_distance(gene1.sequence, gene2.sequence, k)
    return matrix


# lets do the same, but this time use the 'smith_waterman' alignment distance function you wrote. 
def calculate_sw_alignment_distance_matrix(genes, matrix):
    for gene1 in genes:
        for gene2 in genes:
            matrix[gene1.accession][gene2.accession] = smith_waterman(gene1.sequence, gene2.sequence)
    return matrix


distance_matrix = init_distance_matrix(ribosome_genes)
kmer_distance_matrix = calculate_kmer_distance_matrix(ribosome_genes, distance_matrix, 8)
#print('\nkmer distance matrix')
#print_distance_matrix(kmer_distance_matrix)


distance_matrix = init_distance_matrix(ribosome_genes)
#sw_alignment_distance_matrix = calculate_sw_alignment_distance_matrix(ribosome_genes, distance_matrix)
#print('\nsmith waterman alignment score distance matrix')
#print_distance_matrix(sw_alignment_distance_matrix)

#heatmap(kmer_distance_matrix, sw_alignment_distance_matrix)



### Comment on why these two heatmaps look inverted relative to one another
### Comment on the extent of similarity between the Smith-Waterman and k-mer distance approaches, and suggest which may be more appropriate for biological sequences, making references to some knowledge you have about genes.  



##############################################
########### SECTION 2: CLUSTERING ############
##############################################
from utilities import initialise_centroids, average_point, assign_points, plot_kmeans, points_equal, euclidean_distance
# From the heatmaps, it should seem like there are a few clusters in the data
# First, lets convert the pairwise distances to 2D coordinates. 
# This is possible using MDS (link to information)
# After we have transformed the distance matrix to 2D coordinates, we can plot it to see if any clusters are evident.

# provided block
kmer_distances_xy = distance_matrix_to_coordinates_MDS(kmer_distance_matrix)
#sw_distances_xy = distance_matrix_to_coordinates_MDS(sw_alignment_distance_matrix)
#mds_scatterplot(distances_xy)


# seems like there is some clustering happening. 
# Lets use some clustering algorithms to define the clusters. 
# in this manner, we can have an objective way to talk about the patterns in the data.


####### TASK 4: K-means & K-medoids #######
# we are going to use K-means and K-medoids to cluster the data.

# lets implement the k-means algorithm (available in tutorial week 5 worksheet)
# we have provided initialise_centroids, assign_points, and kmeans. 
# write a function which calculates new centroid locations (using the mean)
def calculate_mean_centroids(data, assignments, k):
    centroids = []
    for cluster in range(k):
        points = [point for point, assignment in zip(data, assignments) if assignment == cluster]
        centroids.append(average_point(points))
        
    return centroids


# place calculate_mean_centroids() in the kmeans function below to complete kmeans
def kmeans(data, k):
    d = len(data[0])
    centroids = initialise_centroids(data, k)
    cluster_assignments = assign_points(centroids, data)
    old_centroids = [(0,) * d] * k  # unlikely to be equal to centroids at start
    
    while not points_equal(centroids, old_centroids):
        old_centroids = centroids
        cluster_assignments = assign_points(centroids, data)
        centroids = calculate_mean_centroids(data, cluster_assignments, k)
        #optional
        #plot_kmeans(data, centroids, cluster_assignments, k)
    
    return centroids, cluster_assignments


centroids, cluster_assignments = kmeans(kmer_distances_xy, 2)
#plot_kmeans(distances_xy, centroids, cluster_assignments, 2)

# lets also implement k-medoids while we're at it. The only difference between k-means and k-medoids is the calculate_mean_centroids() step, which will instead be calculate_median_centroids()
# the median can be taken here as the point in the cluster which has smallest cumulative distance to the other points in the cluster
# you can use the provided euclidean_distance() function to calculate distances between points
# write a function which calculates new centroid locations (using the median)
def calculate_median_centroids(data, assignments, k):
    centroids = []
    for cluster in range(k):
        points = [point for point, assignment in zip(data, assignments) if assignment == cluster]

        cumulative_distances = []
        for p1 in points:
            distances = sum([euclidean_distance(p1, p2) for p2 in points])
            cumulative_distances.append(distances)        

        point_index = cumulative_distances.index(min(cumulative_distances))
        centroids.append(points[point_index])
        
    return centroids


# place calculate_mean_centroids() in the kmedoids function below to complete kmedoids
def kmedoids(data, k):
    d = len(data[0])
    centroids = initialise_centroids(data, k)
    cluster_assignments = assign_points(centroids, data)
    old_centroids = [(0,) * d] * k  # unlikely to be equal to centroids at start
    
    while not points_equal(centroids, old_centroids):
        old_centroids = centroids
        cluster_assignments = assign_points(centroids, data)
        centroids = calculate_median_centroids(data, cluster_assignments, k)
        #optional
        #plot_kmeans(data, centroids, cluster_assignments, k)
    
    return centroids, cluster_assignments


centroids, cluster_assignments = kmedoids(kmer_distances_xy, 2)
#plot_kmeans(distances_xy, centroids, cluster_assignments, 2)


# Comment on which algorithm, k-means or k-medoids, performed better on the data.
# Why do you think this is the case? 
# what do you think is the optimal number of clusters?



####### TASK 5: silhouette analysis #######
# It would be nice to have a way to objectively compare clustering performance
# we could compare the selection of k (num clusters) to find the most optimal, or 
# could compare different clustering algorithms to see which performs better. 


# silhouette analysis is one way to do this. 
# Silhouette analysis reports information about how far a datapoint is from the clusters is isn't a part of.
# It produces a plot where we can see the separation between clusters (higher is better) 
# Silhouette coefficients are produced in the range [-1, 1], where 1 indicates the the point is far from other clusters, 0 represents a point which is near the decision boundary between two clusters, and negative values indicate that the point may have been assigned to the wrong cluster. 
# more information is available at https://scikit-learn.org/stable/auto_examples/cluster/plot_kmeans_silhouette_analysis.html


# The cell below gives a visual representation of the silhouette scores for points in each cluster. 
centroids, cluster_assignments = kmeans(kmer_distances_xy, 3)
#plot_kmeans(kmer_distances_xy, centroids, cluster_assignments, 3)
#silhouette_plot(kmer_distances_xy, cluster_assignments, 'kmeans')

#centroids, cluster_assignments = kmedoids(kmer_distances_xy, 3)
#plot_kmeans(kmer_distances_xy, centroids, cluster_assignments, 3)
#silhouette_plot(kmer_distances_xy, cluster_assignments, 'kmedoids')


# the average silhouette score for points in a cluster will be used as a metric to compare 
# clustering performance (higher is better)
# the cell below will run k-means and k-mediods 10 times each and record average silhouette scores for each cluster. 
# This will provide a range of starting centroid locations. 


def silhouette_score_run(distance_matrix, k):
    kmeans_scores = []
    kmedoids_scores = []

    while len(kmeans_scores) < 10:
        try:
            centroids, cluster_assignments = kmeans(distance_matrix, k)
            kmeans_scores.append(silhouette_score(distance_matrix, cluster_assignments))
        except ValueError:
            pass

    while len(kmedoids_scores) < 10:
        try:
            centroids, cluster_assignments = kmedoids(distance_matrix, k)
            kmedoids_scores.append(silhouette_score(distance_matrix, cluster_assignments))
        except ValueError:
            pass

    print('{:^15}{:^15}'.format('kmeans', 'kmedoids'))
    for score1, score2 in zip(kmeans_scores, kmedoids_scores):
        print('{:^15.2f}{:^15.2f}'.format(score1, score2))


print('\nsilhouette scores: k = 2')
silhouette_score_run(kmer_distances_xy, 2)
print('\nsilhouette scores: k = 3')
silhouette_score_run(kmer_distances_xy, 3)



# in light of the silhouette analysis results produced above, which algorithm (k-means/k-medoids) performed better for k=2? 
# comment on why the scores are changing. 



####### TASK 6: divisive clustering algorithm #######

# write a divisive clustering method for the pairwise distance matrix


# Lets assess the performance of your algorithm with the kmeans/kmedoids approaches. 
# Use the silhouette analysis function to get a silhouette score. 
# You can call silhouette_score(distance_matrix, cluster_assignments), where distance_matrix is a pairwise distance matrix, and cluster_assignments are your assignments for each data point. 
# repeat this scoring a number of times if your algorithm is not guaranteed an optimal solution for a given run.



# comment on the performance of the three algorithms and which one was best overall.











