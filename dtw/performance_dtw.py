import random
import time
import openpyxl
from utils import *
from bandDTW import bandDTW
from DTW import dtw
from sparseDTW import sparseDTW
from DC_sparseDTW import DC_sparseDTW
import os

def _sparseDTW(seq_a, seq_b, l1, l2, window_size, res, compare_name):
    '''
    Because sparseDTW may not be the best path, so we look them directly.
    The Performance to be tested:
    1. The final distance
    2. used matrix. It will calculate the cell used to store the distance
    3. elapsed time. The consuming time in each DTW algorithm.
    '''

    file = r"../test/record1.xlsx"
    workbook = openpyxl.load_workbook(file)
    sheet = workbook.active
    data = []

    #Native DTW
    t1 = time.time()
    ndtw = dtw(seq_a, seq_b)
    print("\nNaive DTW")
    print(ndtw())
    t2 = time.time()
    elapsed_time = t2-t1
    distance = ndtw.distance()
    cell_count = ndtw.count_matrix()

    print("Distance:",distance)
    print("Elapsed time:", elapsed_time)
    print("Used matrix:", cell_count)
    data.append([compare_name, l1, l2, elapsed_time, cell_count, distance])

    # BandDTW
    t1 = time.time()
    bdtw = bandDTW(seq_a, seq_b, window_size)
    print("\nBand DTW")
    print(bdtw())
    t2 = time.time()
    elapsed_time = t2 - t1
    distance = bdtw.distance()
    cell_count = bdtw.count_matrix()

    print("Distance:", distance)
    print("Elapsed time:", elapsed_time)
    print("Used matrix:", cell_count)
    data.append([compare_name, l1, l2, elapsed_time, cell_count, distance])

    #sparseDTW
    t1 = time.time()
    sdtw = sparseDTW(seq_a, seq_b, window_size, res)
    print("\nSparse DTW")
    print(sdtw())
    t2 = time.time()
    elapsed_time = t2 - t1
    distance = sdtw.distance()
    cell_count = sdtw.count_matrix()

    print("Distance:", distance)
    print("Elapsed time:", elapsed_time)
    print("Used matrix:", cell_count)
    data.append([compare_name, l1, l2, elapsed_time, cell_count, distance])

    #DC_sparseDTW
    t1 = time.time()
    fdtw = DC_sparseDTW(seq_a, seq_b, window_size, res, 2)
    print("\nDC sparse DTW")
    print(fdtw())
    t2 = time.time()
    elapsed_time = t2 - t1
    distance = fdtw.distance()
    cell_count = fdtw.count_matrix()

    print("Distance:", distance)
    print("Elapsed time:", elapsed_time)
    print("Used matrix:", cell_count)
    data.append([compare_name, l1, l2, elapsed_time, cell_count, distance])
    data.append([])

    for row in data:
        sheet.append(row)
    workbook.save(file)


def per_DC_sparseDTW(seq_a, seq_b, l1, l2, window_size, res, compare_name, radius):
    '''
    Because sparseDTW may not be the best path, so we look them directly.
    The Performance to be tested:
    1. The final distance
    2. used matrix. It will calculate the cell used to store the distance
    3. elapsed time. The consuming time in each DTW algorithm.
    '''

    file = r"../test/record2.xlsx"
    workbook = openpyxl.load_workbook(file)
    sheet = workbook.active
    data = []
    # DC_sparseDTW
    t1 = time.time()
    fdtw = DC_sparseDTW(seq_a, seq_b, window_size, res, radius)
    print("\nDC sparse DTW")
    print(len(seq_a), res, radius)
    fdtw()
    t2 = time.time()
    elapsed_time = t2 - t1
    distance = fdtw.distance()
    cell_count = fdtw.count_matrix()

    print("Distance:", distance)
    print("Elapsed time:", elapsed_time)
    print("Used matrix:", cell_count)
    data.append([compare_name, l1, l2, elapsed_time, cell_count, distance, res, radius])
    data.append([])

    for row in data:
        sheet.append(row)
    workbook.save(file)

def compare_dtw(file1, file2, base_path):

    for x in range(100,2200,200):
        df1 = mat_to_csv(os.path.join(base_path, file1))
        df1 = df1.loc[:,1:x]
        df1 = df1.values
        df2 = mat_to_csv(os.path.join(base_path, file2))
        df2 = df2.loc[:,1:x]
        df2 = df2.values
        _sparseDTW(df1, df2, x, x, x/10, .5, file1.split(".")[0]+ "-" + file2.split(".")[0])

def compare_sparse_dtw(file1, file2, base_path):

    for x in range(100,2200,200):
        df1 = mat_to_csv(os.path.join(base_path, file1))
        df1 = df1.loc[:,1:x]
        df1 = df1.values
        df2 = mat_to_csv(os.path.join(base_path, file2))
        df2 = df2.loc[:,1:x]
        df2 = df2.values
        per_DC_sparseDTW(df1, df2, x, x, x / 10, .3, file1.split(".")[0] + "-" + file2.split(".")[0], 2)
        per_DC_sparseDTW(df1, df2, x, x, x/10, .5, file1.split(".")[0]+ "-" + file2.split(".")[0], 2)
        per_DC_sparseDTW(df1, df2, x, x, x/10, .7, file1.split(".")[0] + "-" + file2.split(".")[0], 2)
        per_DC_sparseDTW(df1, df2, x, x, x / 10, .3, file1.split(".")[0] + "-" + file2.split(".")[0], 5)
        per_DC_sparseDTW(df1, df2, x, x, x / 10, .5, file1.split(".")[0] + "-" + file2.split(".")[0], 5)
        per_DC_sparseDTW(df1, df2, x, x, x / 10, .7, file1.split(".")[0] + "-" + file2.split(".")[0], 5)
        per_DC_sparseDTW(df1, df2, x, x, x / 10, .3, file1.split(".")[0] + "-" + file2.split(".")[0], 10)
        per_DC_sparseDTW(df1, df2, x, x, x / 10, .5, file1.split(".")[0] + "-" + file2.split(".")[0], 10)
        per_DC_sparseDTW(df1, df2, x, x, x / 10, .7, file1.split(".")[0] + "-" + file2.split(".")[0], 10)

if __name__ == '__main__':
    min_num = 10
    max_num = 3581

    # Generate two random file number0s
    for i in range(100):
        file_num_1 = random.randint(min_num, max_num)
        file_num_2 = random.randint(min_num, max_num)

        # Make sure the two file numbers are different
        while file_num_2 == file_num_1:
            file_num_2 = random.randint(min_num, max_num)

        # Create the file names
        try:
            filename1 = 'Q{:04d}.mat'.format(file_num_1)
            filename2 = 'Q{:04d}.mat'.format(file_num_2)
        except Exception:
            continue
        compare_dtw(filename1, filename2, base_path)
     

    df = pd.read_excel('../test/record1.xlsx')
    df = df.iloc[:, 0].drop_duplicates()
    # Get the unique values from the first column
    unique_files = (df.tolist())
    for file_name in unique_files:
        # Split the string by the delimiter
        file_parts = file_name.split('-')
        # Construct the file name by joining the parts with ".mat"
        filename1 = f"{file_parts[0]}.mat"
        filename2 = f"{file_parts[1]}.mat"
        compare_sparse_dtw(filename1, filename2, base_path)




