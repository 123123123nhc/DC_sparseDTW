# DC sparse DTW
DC sparse DTW is an optimized version of Dynamic Time Warping (DTW) that reduces the spatial complexity of the original algorithm while improving its performance for multi-dimensional data comparison. DTW is a widely used technique for finding similarities between time series in various fields, including biological information, audio analysis, and image recognition. However, its quadratic time complexity and space complexity can limit its practical applications. Sparse DTW was a pruning optimization for DTW that reduced its spatial complexity effectively but was still limited to one-dimensional data. DC sparse DTW builds upon sparse DTW by optimizing it using the extended dimension, setting constraints, and implementing divide-and-conquer strategies, resulting in a low spatial complexity algorithm suitable for multi-dimensional data comparison.

# Requirements
#### Python 3.6 or higher
#### numpy==1.16.2
#### scipy==1.3.0
#### pandas==1.1.5
#### matplotlib==3.3.4
#### openpyxl==3.0.7
#### random (part of the Python standard library)
#### math (part of the Python standard library)
#### time (part of the Python standard library)
#### os (part of the Python standard library)

# Installation
To install DC Sparse DTW, you can download or clone the repository and install the required dependencies using github:

```
git clone https://github.com/123123123nhc/DC_sparseDTW.git
```

# Usage
You can use the library by importing it into your Python code:

```
from DC_sparseDTW import DC_sparseDTW

# Generate two random time series
seq1_dim = random.randint(1, 10)
seq1_length_a = random.randint(10, 20)
seq1_length_b = random.randint(10, 20)
seq1_a = [[random.uniform(1, 10) for _ in range(seq1_length_a)] for _ in range(seq1_dim)]
seq1_b = [[random.uniform(1, 10) for _ in range(seq1_length_b)] for _ in range(seq1_dim)]

# Compute the DC sparse DTW distance between the time series
fdtw = DC_sparseDTW(seq1_a, seq1_b, 15, 1)
distance = fdtw()
count_matrix = fdtw.count_matrix()
dtw_path = fdtw.to_array()

print("DC sparse DTW distance:", distance)
print("Count matrix:", count_matrix)
print("DTW path:", dtw_path)
```

Note: You may need to adjust the import statement and other details in the example usage depending on the specific requirements of your project.

# Author
Haocheng Ni
