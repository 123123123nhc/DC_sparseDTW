'''
Author: Haocheng Ni
Date: 2023/3/22
'''

from math import isinf
import random

from numpy import inf
from bandDTW import bandDTW

class sparseDTW(bandDTW):
    '''
    compute the dtw matrix with sparseDTW algorithm. It can save lots of time complexity compared with
    normal dtw because lots of cells in the matrix are blocked (not computed).

    For distance in the matrix, 0 means that the cell between two characters are blocked,
    -1 means that the characters in the cell are the same.
    '''

    def __init__(self, seq_a, seq_b, window_size=inf, res = .5):
        '''
        :param res: the factor to influence the scale of the bin.
        :return: sparse dtw matrix, route path, the score of matching performance.
        '''
        bandDTW.__init__(self, seq_a, seq_b, window_size)
        self.res = res

    def upper_neighbor(self, row, col): #get the coordinates of upper neighbors
        u_neighbor = []
        if row + 1 < self.length_a and abs(row+1-col) < self.w:
            u_neighbor.append((row + 1, col))
        if col + 1 < self.length_b and abs(row-col-1) < self.w:
            u_neighbor.append((row, col + 1))
        if row + 1 < self.length_a and col + 1 < self.length_b:
            u_neighbor.append((row + 1, col + 1))
        return u_neighbor

    def unblocking(self, row, col): #unblock the cell if its upper neighbors are all blocked.
        for upper_neighbor in self.upper_neighbor(row, col):
            if self.matrix[upper_neighbor] != 0:
                return
        for upper_neighbor in self.upper_neighbor(row, col):
            self.set_distance(upper_neighbor[0], upper_neighbor[1])
            self.count += 1

    def get_index(self, sequences, lower_bound, upper_bound): # return the index in each bin
        index_box = []
        seq_length = len(sequences[0])
        for index in range(seq_length):
            index_k = [k for k in range(self.dim) if lower_bound <= sequences[k][index] <= upper_bound]
            if len(index_k) >= self.dim * (1 - self.res):
                index_box.append(index)
        return index_box

    def set_distance(self, row, col):
        if [self.seq_a[k][row] for k in range(self.dim)] != [self.seq_b[k][col] for k in range(self.dim)]:
            super().set_distance(row, col)
        else:
            self.matrix[row, col] = -1

    def set_matrix(self):
        '''
        The first step to generate the sparse dtw matrix. If cells in the matrix are not blocked, it will be changed to
        the euclidean distance of two elements. If two elements are the same, set -1 in the cell.
        :return: the original matrix.
        '''

        if self.w <= abs(self.length_b - self.length_a) :
            raise Exception("Window size is too small. Window size should bigger than %s" % abs(self.length_a - self.length_b))

        lower_bound, upper_bound = 0, self.res

        quantized_seq_a = self.quantize(self.seq_a)
        quantized_seq_b = self.quantize(self.seq_b)
        if self.dim != 1:
            self.seq_a = self.normalize(self.seq_a)
            self.seq_b = self.normalize(self.seq_b)
        self.set_distance(0, 0)
        while lower_bound <= 1 - self.res/2:
            index_set_a = self.get_index(quantized_seq_a, lower_bound, upper_bound)
            index_set_b = self.get_index(quantized_seq_b, lower_bound, upper_bound)
            lower_bound += self.res/2
            upper_bound += self.res
            for index_a in index_set_a:
                for index_b in index_set_b:
                    if not (isinf(self.w) or (max(0, index_a - self.w) <= index_b <= min(self.length_b, index_a + self.w))):
                        continue
                    blocked = False
                    for neighbor in self.lower_neighbor(index_a, index_b):
                        if self.matrix[neighbor] != 0:
                            blocked = False
                            break
                        blocked = True
                    if not blocked:
                        if self.matrix[index_a, index_b] == 0:
                            self.set_distance(index_a, index_b)
                            self.count += 1
        self.set_distance(self.length_a-1, self.length_b-1)

    def update_matrix(self):
        '''
        The second step of generating final matrix. The value in the matrix will be updated to the distance from cell(0,0).
        The cells will be skipped if they are blocked. The distance is 0 if the value in the cell is -1.
        :return: the updated matrix according to sparse dtw algorithm.
        '''
        for row in range(self.length_a):
            for col in range(self.length_b):
                if self.matrix[row,col] == 0:
                    continue
                self.unblocking(row, col)
                if self.matrix[row,col] == -1 and self.lower_neighbor(row, col) != []:
                    self.matrix[row, col] += self.min_cost(row, col) + 1
                else:
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
        if self.matrix[self.min_path(row, col)] == -1:
            return 0
        else:
            return super().min_cost(row, col)

    def find_path(self):
        '''
        get the alignment path by backtracking the smallest cells in the matrix.
        :return: the optimal path. (There may be multiple paths, it will return only one if the optimal path)
        '''
        warping_path = []
        hop = (self.length_a - 1, self.length_b - 1)
        warping_path.append(hop)
        while hop != (0, 0):
            hop = self.min_path(hop[0], hop[1])
            if hop in warping_path:
                self.matrix[warping_path[-1]] = 0
                hop = warping_path[-2]
                warping_path = warping_path[:-2]

            warping_path.append(hop)
        warping_path.reverse()
        return warping_path

class window_sparseDTW(sparseDTW):
    def __init__(self, seq_a, seq_b, window_size, res = .5, window=None):
        sparseDTW.__init__(self, seq_a, seq_b, window_size, res)
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
            if self.w <= abs(self.length_b - self.length_a) :
                raise Exception("Window size is too small. Window size should bigger than %s" % abs(self.length_a - self.length_b))
            lower_bound, upper_bound = 0, self.res
            quantized_seq_a = self.quantize(self.seq_a)
            quantized_seq_b = self.quantize(self.seq_b)
            self.set_distance(0, 0)
            while lower_bound <= 1 - self.res/2:
                index_set_a = self.get_index(quantized_seq_a, lower_bound, upper_bound)
                index_set_b = self.get_index(quantized_seq_b, lower_bound, upper_bound)
                lower_bound += self.res/2
                upper_bound += self.res
                for index_a in index_set_a:
                    for index_b in index_set_b:
                        if not (isinf(self.w) or (max(0, index_a - self.w) <= index_b <= min(self.length_b, index_a + self.w))):
                            continue
                        blocked = False
                        for neighbor in self.lower_neighbor(index_a, index_b):
                            if self.matrix[neighbor] != 0:
                                blocked = False
                                break
                            blocked = True
                        if not blocked:
                            if self.matrix[index_a, index_b] == 0:
                                if (index_a, index_b) in self.window:
                                    self.set_distance(index_a, index_b)
                                    self.count += 1
            self.set_distance(self.length_a-1, self.length_b-1)

if __name__ == '__main__':
    seq1_dim = random.randint(1, 10)
    seq1_a = [[random.randint(1, 10) for _ in range(10)] for _ in range(seq1_dim)]
    seq1_b = [[random.randint(1, 10) for _ in range(10)] for _ in range(seq1_dim)]
    sdtw = sparseDTW(seq1_a, seq1_b,
                     5)
    print(sdtw())
    print(sdtw.lower_neighbor(8,3))
    print(sdtw.to_array())









