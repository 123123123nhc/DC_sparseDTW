import random
from numpy import inf
from bandDTW import bandDTW, window_bandDTW
import numpy as np

class DC_DTW(bandDTW):
    '''
    The dtw algorithm uses divide and conquer algorithm to improve the basic dtw method.
    The method uses window dtw instead of dtw because it used the window to limit the path.

    '''

    def __init__(self, seq_a, seq_b, window_size = inf, radius = 1):
        bandDTW.__init__(self, seq_a, seq_b, window_size)
        self.radius = radius
        self.count = 0
        if self.dim != 1:
            self.seq_a = self.normalize(self.seq_a)
            self.seq_b = self.normalize(self.seq_b)

    def fastDTW(self, seq_a, seq_b, radius):
        '''
        Divide and conquer algorithm. Only calculate path with normal dtw when one of the sequence is short.
        Uses previous path to calculate the window in the subsequence path.
        :param seq_a:
        :param seq_b:
        :param radius: expand the window by path. If the radius is high, it means the window will be expanded broadly.
        Radius should be 1 or 2 in most condition.
        :return: The alignment path
        '''
        if len(seq_a[0]) < radius + 2 or len(seq_b[0]) < radius + 2:
            return window_bandDTW(seq_a, seq_b, self.w)()
        a_merged = self.merge(seq_a)
        b_merged = self.merge(seq_b)
        path = self.fastDTW(a_merged, b_merged, radius)
        window = self.expand_window(path, len(seq_a[0]), len(seq_b[0]), self.radius)
        win_dtw = window_bandDTW(seq_a, seq_b, self.w, window)
        new_path = win_dtw()
        self.matrix = win_dtw.matrix
        self.count = win_dtw.count_matrix()
        return new_path

    def merge(self, seq):
        seq = [[(seq[dim][i] + seq[dim][i+1])/2 for i in range(0, len(seq[0]) - len(seq[0]) % 2, 2)] for dim in range(self.dim)]
        return seq

    def expand_window(self, path, len_x, len_y, radius):
        path_ = set(path)
        for i, j in path:
            for a, b in ((i + a, j + b)
                         for a in range(-radius, radius + 1)
                         for b in range(-radius, radius + 1)):
                path_.add((a, b))

        window_ = set()
        for i, j in path_:
            for a, b in ((i * 2, j * 2), (i * 2, j * 2 + 1),
                         (i * 2 + 1, j * 2), (i * 2 + 1, j * 2 + 1)):
                window_.add((a, b))

        window = []
        start_j = 0
        for i in range(0, len_x):
            new_start_j = None
            for j in range(start_j, len_y):
                if (i, j) in window_:
                    window.append((i, j))
                    if new_start_j is None:
                        new_start_j = j
                elif new_start_j is not None:
                    break
            start_j = new_start_j

        return window

    def distance(self):
        return self.matrix[-1,-1]

    def __call__(self, *args, **kwargs):
        return self.fastDTW(self.seq_a, self.seq_b, self.radius)

if __name__ == '__main__':
    seq1_dim = random.randint(1, 10)
    seq1_a = [[random.uniform(1, 10) for _ in range(10)] for _ in range(seq1_dim)]
    seq1_b = [[random.uniform(1, 10) for _ in range(10)] for _ in range(seq1_dim)]
    fdtw = DC_DTW(seq1_a, seq1_b)
    print(fdtw())
    print(fdtw.count_matrix())
    print(fdtw.to_array())