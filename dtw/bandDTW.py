import random
from numpy import inf
from DTW import dtw
from math import isinf

class bandDTW(dtw):
    '''
    bandDTW : add global constraints on DTW to increase its speed by limiting how far the warping path
    may stray from the diagonal of the warping matrix
    '''
    def __init__(self, seq_a, seq_b, window_size=inf):
        dtw.__init__(self, seq_a, seq_b)
        self.w = window_size

    def set_matrix(self):
        '''
        The first step to generate the sparse dtw matrix. If cells in the matrix are not blocked, it will be changed to
        the euclidean distance of two elements. If two elements are the same, set -1 in the cell.
        :return: the original matrix.
        '''
        if self.w <= abs(self.length_b - self.length_a) :
            raise Exception("Window size is too small. Window size should bigger than %s" % abs(self.length_a - self.length_b))
        if self.dim != 1:
            self.seq_a = self.normalize(self.seq_a)
            self.seq_b = self.normalize(self.seq_b)
        for index_a in range(self.length_a):
            for index_b in range(self.length_b):
                if (isinf(self.w) or (max(0, index_a - self.w) <= index_b <= min(self.length_b, index_a + self.w))):
                    self.set_distance(index_a, index_b)
                    self.count += 1

    def lower_neighbor(self, row, col): #get the coordinates of lower neighbors
        l_neighbor = []
        if row - 1 >= 0 and abs(row-1-col) < self.w:
            l_neighbor.append((row-1, col))
        if col - 1 >= 0 and abs(row-col+1) < self.w:
            l_neighbor.append((row, col-1))
        if row - 1 >= 0 and col - 1 >= 0:
            l_neighbor.append((row-1, col-1))
        return l_neighbor


class window_bandDTW(bandDTW):
    def __init__(self, seq_a, seq_b, window_size, window = None):
        bandDTW.__init__(self,seq_a, seq_b, window_size)
        self.window = window

    def set_matrix(self):
        '''
        The first step to generate the sparse dtw matrix. If cells in the matrix are not blocked, it will be changed to
        the euclidean distance of two elements. If two elements are the same, set -1 in the cell.
        :return: the original matrix.
        '''
        if self.window is None:
            super().set_matrix()
        else:
            if self.w <= abs(self.length_b - self.length_a):
                raise Exception(
                    "Window size is too small. Window size should bigger than %s" % abs(self.length_a - self.length_b))
            if self.dim != 1:
                self.seq_a = self.normalize(self.seq_a)
                self.seq_b = self.normalize(self.seq_b)
            for index_a in range(self.length_a):
                for index_b in range(self.length_b):
                    if (isinf(self.w) or (max(0, index_a - self.w) <= index_b <= min(self.length_b, index_a + self.w))):
                        if (index_a, index_b) in self.window:
                            self.set_distance(index_a, index_b)
                            self.count += 1

    def lower_neighbor(self, row, col):
        l_neighbors = super().lower_neighbor(row, col)
        if self.window is None:
            return l_neighbors
        return [l_neighbor for l_neighbor in l_neighbors if l_neighbor in self.window]


if __name__ == '__main__':
    seq1_dim = random.randint(1, 10)
    seq1_a = [[random.uniform(1, 10) for _ in range(10)] for _ in range(seq1_dim)]
    seq1_b = [[random.uniform(1, 10) for _ in range(10)] for _ in range(seq1_dim)]
    bdtw = bandDTW(seq1_a, seq1_b, 5)
    print(bdtw())
    print(bdtw.to_array())


