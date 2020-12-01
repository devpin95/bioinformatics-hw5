import argparse
from kmeans import K_Means

parser = argparse.ArgumentParser()
parser.add_argument("exp_data", help="name of a file containing expression data")
parser.add_argument("k", help="an integer indicating k, the number of clusters to form", type=int)
parser.add_argument("--iter", help="the number of iterations to train", type=int)
args = parser.parse_args()

exp_data = open(args.exp_data)
k = args.k
iter = args.iter if args.iter is not None else 5

# print("k = " + str(k) + " on file " + args.exp_data)

genes = []
for line in exp_data:
    gene = {}
    gene_str = line.split('\t')

    gene['id'] = gene_str[0]
    gene['description'] = gene_str[1]

    for i in range(0, len(gene_str[2:])):
        if gene_str[i] == '':
            gene_str[i] = '0.0'

    gene['exp'] = [float(val) for val in gene_str[2:] if val != '']

    genes.append(gene)

knn = K_Means(k, genes, iter)
knn.train()
knn.describe_clusters()
