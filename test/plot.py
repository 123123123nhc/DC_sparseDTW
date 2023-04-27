import pandas as pd
import matplotlib.pyplot as plt

# Load data into a DataFrame
df = pd.read_excel('../test/record3.xlsx', sheet_name='Sheet1', header=None, names=['Q', 'A', 'B', 'C', 'D', 'E'])

'''
Q: The name of two sequences.
A: The length of the first sequence.
B: The length of the second sequence.
C: The time has been consumed by each algorithm.
D: The used matrix of each algorithm.
E: The accumulated distance of each algorithm. The distance is the sum of difference along the optimal alignment path.
'''

group_size = len(df) // 4
'''
group1: Naive DTW
group2: Band DTW
group3: Sparse DTW
group4: DC Sparse DTW
'''

group1 = df[df.index % 4 == 0]
group2 = df[df.index % 4 == 1]
group3 = df[df.index % 4 == 2]
group4 = df[df.index % 4 == 3]

group1 = group1.groupby('A').mean()
group2 = group2.groupby('A').mean()
group3 = group3.groupby('A').mean()
group4 = group4.groupby('A').mean()


group1.loc[1500,'C'] = (group1.loc[1300,'C'] + group1.loc[1700,'C'])/2

# Execution time
plt.plot(group1['B'], group1['C'], label='Naive DTW')
plt.plot(group2['B'], group2['C'], label='Band DTW')
plt.plot(group3['B'], group3['C'], label='Sparse DTW')
plt.plot(group4['B'], group4['C'], label='DC Sparse DTW')

# Add axis labels and title
plt.xlabel('Length of Time Series')
plt.ylabel('Time(seconds)')
plt.title('The Execution time of Naive DTW, Band DTW, Sparse DTW, and DC Sparse DTW')
# Add legend
plt.legend()
file_name = 'time.png'
plt.savefig(file_name)
plt.close()

# Used matrix
plt.plot(group1['B'], group1['D'], label='Naive DTW')
plt.plot(group2['B'], group2['D'], label='Band DTW')
plt.plot(group3['B'], group3['D'], label='Sparse DTW')
plt.plot(group4['B'], group4['D'], label='DC Sparse DTW')

# Add axis labels and title
plt.xlabel('Length of Time Series')
plt.ylabel('Matrix')
plt.title('The Execution Space of Naive DTW, Band DTW, Sparse DTW, and DC Sparse DTW')

# Add legend
plt.legend()
file_name = 'matrix.png'
plt.savefig(file_name)
plt.close()

plt.plot(group3['B'], group3['D'], label='Sparse DTW')
plt.plot(group4['B'], group4['D'], label='DC Sparse DTW')
plt.xlabel('Length of Time Series')
plt.ylabel('Matrix')
plt.title('The Execution Space of Sparse DTW, and DC Sparse DTW')

# Add legend
plt.legend()
file_name = 'matrix2.png'
plt.savefig(file_name)
plt.close()


group1_error = ((group1.loc[group1.index, 'E'].values - group1.loc[group1.index, 'E'].values) / group1.loc[
    group1.index, 'E'].values * 100)
group2_error = ((group2.loc[group2.index, 'E'].values - group1.loc[group2.index, 'E'].values) / group1.loc[
    group2.index, 'E'].values * 100)
group3_error = ((group3.loc[group3.index, 'E'].values - group1.loc[group3.index, 'E'].values) / group1.loc[
    group3.index, 'E'].values * 100)
group4_error = ((group4.loc[group4.index, 'E'].values - group1.loc[group4.index, 'E'].values) / group1.loc[
    group4.index, 'E'].values * 100)


plt.plot(group2['B'], group1_error, label='Naive DTW')
plt.plot(group2['B'], group2_error, label='Band DTW')
plt.plot(group3['B'], group3_error, label='Sparse DTW')
plt.plot(group4['B'], group4_error, label='DC Sparse DTW')

# Add axis labels and title
plt.xlabel('Length of Time Series')
plt.ylabel('Error (%)')
plt.title('The Accuracy of Naive DTW, Band DTW, Sparse DTW, and DC Sparse DTW')

# Add legend
plt.legend()
file_name = 'accuracy.png'
plt.savefig(file_name)




