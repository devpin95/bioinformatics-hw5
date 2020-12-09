import math
from centron import Centron

NUM_TRAINING_GENES = 3

class K_Means:
    def __init__(self, k, genes, iter):
        self.k = k  # number of clusters
        self.genes = genes  # a list of gene data points
        self.num_genes = len(genes)  # the number of genes we are training on
        self.centrons = []  # list of centrons to train
        self.dim = len(genes[0]['exp'])  # the length of the vector describing each gene
        self.iter = iter  # the number of iterations fort raining

        self.init_centrons()

    def init_centrons(self):
        # initialize the list of centrons that we are going to train
        # each centron gets 3 genes to initialize it's point
        # the first centron gets the first 3 genes, the next centron gets the next 3 genes, and so on

        genes_per_centron = math.floor(self.num_genes / self.k)

        # loop through each centron
        for i in range(0, self.k):
            centron = Centron(self.dim)  # initialize a centron with a self.dim length vector of zeros
            self.centrons.append(centron)  # add the new centron to the list

        # loop through the genes to initialize the centrons
        for i in range(0, self.num_genes):
            # use the number of training genes to determine which centron to add the current gene to
            centron = math.floor(i/NUM_TRAINING_GENES)

            # if we are looking at a centron outside the list of our list of centrons, break
            if centron >= len(self.centrons):
                break

            # otherwise, add the gene to the neighborhood of the current centron
            self.centrons[centron].add_neighbor(self.genes[i]['id'])

        # make the centrons update their position based on their neighborhood
        self.update_centrons()

    def train(self):
        # main training loop for clustering

        # loop the number of times set by self.iter
        for epoch in range(0, self.iter):
            # print("EPOCH " + str(epoch) + '-------------------------------------------------------------------------')
            self.reset_centrons()  # reset each centron's neighborhood
            self.find_nearest_centron()  # go through the points and find the centron's new neighborhood
            self.update_centrons()  # update the centron again with it's new neighbors

    def update_centrons(self):
        # update each centron based on the genes in it's neighborhood

        # loop through the centrons
        for i in range(0, len(self.centrons)):
            # get a list of the genes closest to the current centron, and the number of neighbors in the list
            neighbors = self.centrons[i].neighbors
            num_neighbors = len(neighbors)

            # break from the loop if the centron has 0 neighbors, we don't want to divide by 0
            if num_neighbors == 0:
                break

            # create a new vector of zeros that will hold the sum of each point in the neighborhood
            new_vector = [0.0] * self.dim

            # loop through the neighbors
            for neighbor in neighbors:
                # get the vector for this gene
                vector = next(d['exp'] for d in self.genes if d['id'] == neighbor)

                # loop through the gene's vector and add it to the new vector, divide by the number of neighbors to avg
                for j in range(0, self.dim):
                    new_vector[j] += vector[j] / num_neighbors

            # print("\tCentron " + str(i) + " moved " + str(self.euclidean_distance(new_vector, self.centrons[i].vector)))

            # update the centron with it's new center of mass
            self.centrons[i].update_point(new_vector)

    def reset_centrons(self):
        # empty each centron's neighbor list so that we can find it's new neighbors after updating the centrons
        for centron in self.centrons:
            centron.clear_neighbors()

    def find_nearest_centron(self):
        # find the centron closest to a gene

        # loop through the gene list
        for gene in self.genes:
            # create a vector of distances to centrons, initialize to inf in case something goes wrong
            dist_vector = [math.inf] * self.k

            # loop through the centrons
            for i in range(0, self.k):
                # find the euclidean distance from the gene to the centron
                dist_vector[i] = self.euclidean_distance(self.centrons[i].vector, gene['exp'])

            # the closest centron to this gene will be the min val in the vector
            closest_centron = dist_vector.index(min(dist_vector))

            # add the gene to the neighborhood of the centron it is closest to
            self.centrons[closest_centron].add_neighbor(gene['id'])

    def euclidean_distance(self, data, centron):
        # finds the distance between two points
        # euclidean distance for this data will be
        # sqrt[ (x1,1-x1,2)^2 + (x2,1-x2,2)^2 + (x3,1-x3,2)^2 + ... + (x79,1-x79,2)^2 ]
        # for xi,j where i in [0, self.dim) and j = 1 is the centron vector and j = 2 is the gene vector

        esum = 0.0  # set the sum to 0, this is the value inside the sqrt
        for i in range(0, self.dim):
            # add the difference of the 2 points and square it
            esum += (centron[i] - data[i]) ** 2

        # return the euclidean distance
        return math.sqrt(esum)

    def describe_clusters(self):
        # print the necessary information about the clusters

        # sort the clusters
        sorted_centrons = self.sort_clusters()

        # loop through the sorted centrons
        for centron in sorted_centrons:
            # sort the genes in the neighborhood of the centron
            sorted_neighbors = self.sort_centron_neighbors(centron['centron_index'])

            # loop through the sorted genes
            for gene_data in sorted_neighbors:
                # get the gene metadata (id and description) and print it
                gene = next(d for d in self.genes if d['id'] == gene_data['gene'])
                print(gene['id'], gene['description'])

            # get the centron's vector data
            cluster_vector = self.centrons[centron['centron_index']].vector

            # loop through the centron vector and print the rounded value
            for val in cluster_vector:
                print(str(round(val, 4)) + '\t', end='')
            print('\n')

    def sort_clusters(self):
        # sort the centrons by their average expression rate

        # list to hold the centrons
        centrons = []

        # loop through the centrons
        for i in range(0, self.k):
            # create a dict with the centron id (it's index) and it's average expression rate
            data = {'centron_index': i, 'avg_exp': self.average_exprate(self.centrons[i].vector)}
            centrons.append(data)

        # return the list sorted by the lambda (which sorts on the avg_exp key in the dict)
        return sorted(centrons, key=lambda i: i['avg_exp'])

    def sort_centron_neighbors(self, centron):
        # sort the genes in the neighborhood of a centron by their average expression rate

        # get the list of neighbors for the centron
        neighbors = self.centrons[centron].neighbors

        # create a list to hold the centrons to sort later
        neighbor_data = []

        # loop through the neighbors
        for neighbor in neighbors:
            # get the vector for the neighbor
            vector = next(d['exp'] for d in self.genes if d['id'] == neighbor)

            # create a dict to hold metadata about the gene (the id and it's average expression rate)
            data = {'gene': neighbor, 'avg_exp': self.average_exprate(vector)}
            neighbor_data.append(data)

        # return the list sorted by the lambda (which sorts on the avg_exp key in the dict)
        return sorted(neighbor_data, key=lambda i: i['avg_exp'])

    def average_exprate(self, vector):
        # finds the average expression rate of a vector
        # sums the entire vector and divides by it's length

        esum = 0
        for val in vector:
            esum += val

        return esum / len(vector)
