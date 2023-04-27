import random
from scipy.sparse import lil_matrix
import numpy as np

'''
It is the root class of all DTW method. Although some of the functions may not be used in this class, it is useful in 
the other methods. dtw is the method without window constraint. window_dtw is the method that can be limited by 
the parameter, window. 
'''

class dtw:
    '''
    The parent class of bandDTW and sparseDTW. The 1D DTW or multiple dimensional dtw can be calculated by this class.
    1D DTW: normal DTW algorithm;
    MD DTW: accumulate the distance of each dimensions by l1 normalization.
    All the functions will be used in bandDTW and sparseDTW.
    '''
    def __init__(self, seq_a, seq_b):
        '''
        :param seq_a: any sequence with a fixed order.
        :param seq_b: Another sequence with a fixed order.
        '''
        self.seq_a = [seq_a] if np.array(seq_a).ndim == 1 else seq_a
        self.seq_b = [seq_b] if np.array(seq_b).ndim == 1 else seq_b
        assert len(self.seq_a) == len(self.seq_b), "The dimension of two sequences are not equal!"
        self.dim = len(self.seq_a)
        self.length_a, self.length_b = len(self.seq_a[0]), len(self.seq_b[0])
        self.matrix = lil_matrix((self.length_a, self.length_b), dtype=np.float64)
        self.count = 0

    def euc_distance(self, char_a, char_b):
        return (char_a - char_b) ** 2

    def quantize(self, sequence):
        return [[(sequence[dim][i] - min(sequence[dim])) / (max(sequence[dim]) - min(sequence[dim])) for i in
                 range(len(sequence[dim]))] for dim in range(self.dim)]

    def normalize(self, sequence):  # return the value of sequence to make sure that the value is within (0,1)
        return [[(sequence[dim][i] - np.mean(sequence[dim])) / np.std(sequence[dim]) for i in
                 range(len(sequence[dim]))] for dim in range(self.dim)]

    def set_distance(self, row, col):
        dis = 0
        for i in range(self.dim):
            dis += self.euc_distance(self.seq_a[i][row], self.seq_b[i][col])
        self.matrix[row,col] = dis

    def lower_neighbor(self, row, col): #get the coordinates of lower neighbors
        l_neighbor = []
        if row - 1 >= 0:
            l_neighbor.append((row-1, col))
        if col - 1 >= 0:
            l_neighbor.append((row, col-1))
        if row - 1 >= 0 and col - 1 >= 0:
            l_neighbor.append((row-1, col-1))
        return l_neighbor


    def set_matrix(self):
        '''
        The first step to generate the sparse dtw matrix. If MD DTW: quantize the sequences first.
        :return: the original matrix.
        '''
        if self.dim != 1:
            self.seq_a = self.normalize(self.seq_a)
            self.seq_b = self.normalize(self.seq_b)
        for index_a in range(self.length_a):
            for index_b in range(self.length_b):
                self.set_distance(index_a, index_b)
                self.count += 1

    def update_matrix(self):
        '''
        The second step of generating final matrix. The value in the matrix will be updated to the distance from cell(0,0).

        '''
        for row in range(self.length_a):
            for col in range(self.length_b):
                self.matrix[row, col] += self.min_cost(row, col)

    def min_path(self, row, col):#find the coordinate of lower neighbors which contains the lowest value.
        unblocked_lower_neighbors = [lower_neighbor for lower_neighbor in self.lower_neighbor(row, col)
                                     if self.matrix[lower_neighbor] != 0]
        if unblocked_lower_neighbors == []:
            return (row, col)
        min_neighbor = unblocked_lower_neighbors[0]
        for neighbor in unblocked_lower_neighbors:
            if self.matrix[neighbor] <= self.matrix[min_neighbor]:
                min_neighbor = neighbor
        return min_neighbor

    def min_cost(self, row, col): #find the value in the min_path coordinate.
        if (row, col) == self.min_path(row, col):
            return 0
        else:
            return self.matrix[self.min_path(row, col)]

    def find_path(self):
        '''
        get the alignment path by backtracking the smallest cells in the matrix.
        :return: the optimal path. (There may be multiple paths, it will return only one if the optimal path)
        '''
        warping_path = []
        hop = (self.length_a - 1, self.length_b - 1)
        warping_path.append(hop)
        while hop != (0, 0):
            hop = self.min_path(hop[0],hop[1])
            warping_path.append(hop)
        warping_path.reverse()
        return warping_path

    def to_array(self): #matrix visualization
        round_matrix = self.matrix
        for row in range(self.length_a):
            for col in range(self.length_b):
                round_matrix[row, col] = np.round(round_matrix[row, col], 2)

        return round_matrix.toarray()

    def __call__(self, *args, **kwargs):
        self.set_matrix()
        self.update_matrix()
        return self.find_path()

    def count_matrix(self):
        return self.count

    def distance(self):
        return self.matrix[-1, -1]

class window_DTW(dtw):
    '''
    It is an additional version of dtw with a new parameter window.
    window: a List of possible root. The road will not pass if the cell not in window. It is used in divide and conquer
    algorithm.
    '''
    def __init__(self, seq_a, seq_b, window = None):
        dtw.__init__(self, seq_a, seq_b)
        self.window = window

    def set_matrix(self):
        '''
        Additional set matrix method considering the window.
        '''
        if self.window is None:
            super().set_matrix()
        else:
            for index_a in range(self.length_a):
                for index_b in range(self.length_b):
                    if (index_a, index_b) in self.window:
                        self.set_distance(index_a, index_b)
                        self.count += 1

    def update_matrix(self):
        '''
        The second step of generating final matrix. The value in the matrix will be updated to the distance from cell(0,0).

        '''
        if self.window is None:
            super().update_matrix()
        else:
            for (row, col) in self.window:
                self.matrix[row, col] += self.min_cost(row, col)

    def lower_neighbor(self, row, col):
        l_neighbors = super().lower_neighbor(row, col)
        if self.window is None:
            return l_neighbors
        return [l_neighbor for l_neighbor in l_neighbors if l_neighbor in self.window]

if __name__ == '__main__':
    seq1_dim = random.randint(1, 10)
    seq1_length_a = random.randint(10, 30)
    seq1_length_b = random.randint(10, 30)
    seq1_a = [random.uniform(1, 10) for _ in range(seq1_length_a)]
    seq1_b = [random.uniform(1, 10) for _ in range(seq1_length_b)]
    dtw = dtw(seq1_a, seq1_b)
    print(dtw.seq_a)
    print("")
    print(dtw())
    print(dtw.to_array())
