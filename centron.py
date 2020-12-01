class Centron:
    def __init__(self, dim):
        self.vector = [0.0] * dim
        self.neighbors = []

    def add_neighbor(self, neighbor):
        if neighbor not in self.neighbors:
            self.neighbors.append(neighbor)

    def update_point(self, point):
        self.vector = point.copy()

    def clear_neighbors(self):
        self.neighbors = []