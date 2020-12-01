import math
from centron import Centron

NUM_TRAINING_GENES = 3

class K_Means:
    def __init__(self, k, genes, iter):
        self.k = k
        self.genes = genes
        self.num_genes = len(genes)
        self.centrons = []
        self.dim = len(genes[0]['exp'])
        self.iter = iter

        self.init_centrons()

    def init_centrons(self):
        genes_per_centron = math.floor(self.num_genes / self.k)

        for i in range(0, self.k):
            centron = Centron(self.dim)
            self.centrons.append(centron)

        for i in range(0, self.num_genes):
            centron = math.floor(i/NUM_TRAINING_GENES)
            if centron >= len(self.centrons):
                break
            self.centrons[centron].add_neighbor(self.genes[i]['id'])

        self.update_centrons()

    def train(self):
        for epoch in range(0, self.iter):
            # print("EPOCH " + str(epoch) + '---------------------------------------------------------------------------')
            self.reset_centrons()
            self.find_nearest_centron()
            self.update_centrons()

    def update_centrons(self):
        for i in range(0, len(self.centrons)):
            neighbors = self.centrons[i].neighbors
            num_neighbors = len(neighbors)

            if num_neighbors == 0:
                break

            new_vector = [0.0] * self.dim
            for neighbor in neighbors:
                vector = next(d['exp'] for d in self.genes if d['id'] == neighbor)
                for j in range(0, self.dim):
                    new_vector[j] += vector[j] / num_neighbors

            # print("\tCentron " + str(i) + " moved " + str(self.euclidean_distance(new_vector, self.centrons[i].vector)))
            self.centrons[i].update_point(new_vector)

    def reset_centrons(self):
        for centron in self.centrons:
            centron.clear_neighbors()

    def find_nearest_centron(self):
        for gene in self.genes:
            dist_vector = [math.inf] * self.k
            for i in range(0, self.k):
                dist_vector[i] = self.euclidean_distance(self.centrons[i].vector, gene['exp'])

            closest_centron = dist_vector.index(min(dist_vector))
            self.centrons[closest_centron].add_neighbor(gene['id'])

    def euclidean_distance(self, data, centron):
        sum = 0.0
        for i in range(0, self.dim):
            sum += (centron[i] - data[i]) ** 2
        return math.sqrt(sum)

    def describe_clusters(self):
        sorted_centrons = self.sort_clusters()

        for centron in sorted_centrons:
            sorted_neighbors = self.sort_centron_neighbors(centron['centron_index'])

            for gene_data in sorted_neighbors:
                gene = next(d for d in self.genes if d['id'] == gene_data['gene'])
                print(gene['id'], gene['description'])

            cluster_vector = self.centrons[centron['centron_index']].vector

            for val in cluster_vector:
                print(str(round(val, 4)) + '\t', end='')
            print('\n')

    def sort_clusters(self):
        centrons = []
        for i in range(0, self.k):
            data = {'centron_index': i, 'avg_exp': self.average_exprate(self.centrons[i].vector)}
            centrons.append(data)

        return sorted(centrons, key=lambda i: i['avg_exp'])

    def sort_centron_neighbors(self, centron):
        neighbors = self.centrons[centron].neighbors
        neighbor_data = []
        for neighbor in neighbors:
            vector = next(d['exp'] for d in self.genes if d['id'] == neighbor)
            data = {'gene': neighbor, 'avg_exp': self.average_exprate(vector)}
            neighbor_data.append(data)

        return sorted(neighbor_data, key=lambda i: i['avg_exp'])

    def average_exprate(self, vector):
        sum = 0
        for val in vector:
            sum += val

        return sum / len(vector)
