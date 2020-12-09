class Centron:
    # holds information about a centron for kmeans clustering
    def __init__(self, dim):
        self.vector = [0.0] * dim  # the length of the vector, initialized to zeros
        self.neighbors = []  # a list of genes that are closest to this centron

    def add_neighbor(self, neighbor):
        # add a gene to the neighborhood of  this centron
        # first, check if the gene is already in the neighborhood before adding it again
        if neighbor not in self.neighbors:
            self.neighbors.append(neighbor)

    def update_point(self, point):
        # updates the centron's vector to the point passed in
        self.vector = point.copy()

    def clear_neighbors(self):
        # empties the neighbor list so that it can be updated with after the new center of mass is found
        self.neighbors = []
