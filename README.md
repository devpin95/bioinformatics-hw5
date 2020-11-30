# Bioinformatics Homework 5 [CSCI 4800]

## Finding Gene Cluster

You are to implement the *k-means clustering* algorithm and apply it to a small gene-expression data set. Your code should read gene expression measurements from a file, do the clustering for a specified value ofk, and then print out a description of each cluster.

Your program should use Euclidean distance when comparing gene expression profiles to cluster means.

To initialize your cluster means, you should use the following procedure. For each of the *k* clusters, you should select three genes from the input file starting with the first gene in the file and going in order. The initial mean for a given cluster should be the mean vector of these three gene expression profiles. For example, ifk = 3, the first three genes in the file are used to initialize the first cluster, the fourth through sixth genes are used to initialize the second cluster, and the seventh through ninth genes are used to initialize the third cluster.

You can choose to write any programming languages. Please make sure to name your program “Cluster”.

## Input Data
Your program should take the following as command-line arguments:

1. the name of a file containing expression data
2. an integer indicating *k*, the number of clusters to form

You should test your program using the following two files:
*  a set of 12 yeast genes
*  a set of 52 yeast genes

You can assume that formats of any files that your program handles will be the same as these. Each line of the file describes a single gene, and tabs are used to delimit the individual columns of a line. The first column lists an identifier for each gene. The second column lists a common name for the gene and a description of its function. The remaining columns list expression values for the gene under various conditions. These values are log expression ratios.

There are some missing values in the input file (i.e., in some cases there is not a number between the tabs that separate the columns). These missing values indicate cases for which the expression measurement process failed for some reason. You should handle these cases by treating them as equivalent to values of 0 (i.e., expression ratios of 1, no change from the baseline condition).

Your program should not use the gene descriptions in doing the clustering. However, it may be informative to see them in the output. These genes belong to four categories for which the genes in each should exhibit fairly similar expression profiles.

## Output
Your program should output a description of each cluster. In particular, for each cluster your program should list the identifiers and descriptions of the genes that belong to it. You should list these genes, one per line, on consecutive lines. Additionally, after listing the genes that belong to a given cluster, you should list the coordinates of the cluster mean on a separate line. There should be a blank line between each of these cluster descriptions.

You should print the clusters __and__ the genes within each cluster __in order__ as follows. The clusters should be ordered by the average expression ratio, across the various measurements, of the genes they contain (from smallest to largest). Similarly, the genes within a cluster should be ordered by their average expression ratio (from smallest to largest).

