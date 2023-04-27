import pandas as pd
import matplotlib.pyplot as plt

# Load data into a DataFrame
df = pd.read_excel('../test/record5.xlsx', sheet_name='Sheet1', header=None, names=['Q', 'A', 'B', 'C', 'D', 'E', 'F', 'G'])
df2 = pd.read_excel('../test/record3.xlsx', sheet_name='Sheet1', header=None, names=['Q', 'A', 'B', 'C', 'D', 'E'])

'''
Q: The name of two sequences.
A: The length of the first sequence.
B: The length of the second sequence.
C: The time has been consumed by each algorithm.
D: The used matrix of each algorithm.
E: The accumulated distance of each algorithm. The distance is the sum of difference along the optimal alignment path.
F: The res of DC_sparse DTW
G: The radius of DC_sparse DTW
'''

group_size = len(df) // 9
'''
group1: res = .3, radius = 2
group2: res = .5, radius = 2
group3: res = .7, radius = 2
group4: res = .3, radius = 5
group5: res = .5, radius = 5
group6: res = .7, radius = 5
group7: res = .3, radius = 10
group8: res = .5, radius = 10
group9: res = .7, radius = 10

'''

group0 = df2[df2.index % 4 == 0]
group1 = df[df.index % 9 == 0]
group2 = df[df.index % 9 == 1]
group3 = df[df.index % 9 == 2]
group4 = df[df.index % 9 == 3]
group5 = df[df.index % 9 == 4]
group6 = df[df.index % 9 == 5]
group7 = df[df.index % 9 == 6]
group8 = df[df.index % 9 == 7]
group9 = df[df.index % 9 == 8]

group0 = group0.groupby('A').mean()
group1 = group1.groupby('A').mean()
group2 = group2.groupby('A').mean()
group3 = group3.groupby('A').mean()
group4 = group4.groupby('A').mean()
group5 = group5.groupby('A').mean()
group6 = group6.groupby('A').mean()
group7 = group7.groupby('A').mean()
group8 = group8.groupby('A').mean()
group9 = group9.groupby('A').mean()

# Execution time
plt.plot(group1['B'], group1['C'], label='DC sparse DTW(res = .3, radius = 2)')
plt.plot(group2['B'], group2['C'], label='DC sparse DTW(res = .5, radius = 2)')
plt.plot(group3['B'], group3['C'], label='DC sparse DTW(res = .8, radius = 2)')
plt.plot(group4['B'], group4['C'], label='DC sparse DTW(res = .3, radius = 5)')
plt.plot(group5['B'], group5['C'], label='DC sparse DTW(res = .5, radius = 5)')
plt.plot(group6['B'], group6['C'], label='DC sparse DTW(res = .8, radius = 5)')
plt.plot(group7['B'], group7['C'], label='DC sparse DTW(res = .3, radius = 10)')
plt.plot(group8['B'], group8['C'], label='DC sparse DTW(res = .5, radius = 10)')
plt.plot(group9['B'], group9['C'], label='DC sparse DTW(res = .8, radius = 20)')

# Add axis labels and title
plt.xlabel('Length of Time Series')
plt.ylabel('Time(seconds)')
plt.title('The Execution time of DC Sparse DTW')
# Add legend
plt.legend()
file_name = 'DC_time.png'
plt.savefig(file_name)
plt.close()

# Used matrix
plt.plot(group1['B'], group1['D'], label='DC sparse DTW(res = .3, radius = 2)')
plt.plot(group2['B'], group2['D'], label='DC sparse DTW(res = .5, radius = 2)')
plt.plot(group3['B'], group3['D'], label='DC sparse DTW(res = .8, radius = 2)')
plt.plot(group4['B'], group4['D'], label='DC sparse DTW(res = .3, radius = 5)')
plt.plot(group5['B'], group5['D'], label='DC sparse DTW(res = .5, radius = 5)')
plt.plot(group6['B'], group6['D'], label='DC sparse DTW(res = .8, radius = 5)')
plt.plot(group7['B'], group7['D'], label='DC sparse DTW(res = .3, radius = 20)')
plt.plot(group8['B'], group8['D'], label='DC sparse DTW(res = .5, radius = 20)')
plt.plot(group9['B'], group9['D'], label='DC sparse DTW(res = .8, radius = 20)')

# Add axis labels and title
plt.xlabel('Length of Time Series')
plt.ylabel('Matrix')
plt.title('The Execution Space of DC Sparse DTW')

# Add legend
plt.legend()
file_name = 'DC_matrix.png'
plt.savefig(file_name)
plt.close()

group1_error = ((group1.loc[group1.index, 'E'].values - group0.loc[group0.index, 'E'].values) / group0.loc[
    group0.index, 'E'].values * 100)
group2_error = ((group2.loc[group1.index, 'E'].values - group0.loc[group0.index, 'E'].values) / group0.loc[
    group0.index, 'E'].values * 100)
group3_error = ((group3.loc[group1.index, 'E'].values - group0.loc[group0.index, 'E'].values) / group0.loc[
    group0.index, 'E'].values * 100)
group4_error = ((group4.loc[group1.index, 'E'].values - group0.loc[group0.index, 'E'].values) / group0.loc[
    group0.index, 'E'].values * 100)
group5_error = ((group5.loc[group1.index, 'E'].values - group0.loc[group0.index, 'E'].values) / group0.loc[
    group0.index, 'E'].values * 100)
group6_error = ((group6.loc[group1.index, 'E'].values - group0.loc[group0.index, 'E'].values) / group0.loc[
    group0.index, 'E'].values * 100)
group7_error = ((group7.loc[group1.index, 'E'].values - group0.loc[group0.index, 'E'].values) / group0.loc[
    group0.index, 'E'].values * 100)
group8_error = ((group8.loc[group1.index, 'E'].values - group0.loc[group0.index, 'E'].values) / group0.loc[
    group0.index, 'E'].values * 100)
group9_error = ((group9.loc[group1.index, 'E'].values - group0.loc[group0.index, 'E'].values) / group0.loc[
    group0.index, 'E'].values * 100)


plt.plot(group1['B'], group1_error, label='DC sparse DTW(res = .3, radius = 2)')
plt.plot(group2['B'], group2_error, label='DC sparse DTW(res = .5, radius = 2)')
plt.plot(group3['B'], group3_error, label='DC sparse DTW(res = .8, radius = 2)')
plt.plot(group4['B'], group4_error, label='DC sparse DTW(res = .3, radius = 5)')
plt.plot(group5['B'], group5_error, label='DC sparse DTW(res = .5, radius = 5)')
plt.plot(group6['B'], group6_error, label='DC sparse DTW(res = .8, radius = 5)')
plt.plot(group7['B'], group7_error, label='DC sparse DTW(res = .3, radius = 20)')
plt.plot(group8['B'], group8_error, label='DC sparse DTW(res = .5, radius = 20)')
plt.plot(group9['B'], group9_error, label='DC sparse DTW(res = .8, radius = 20)')

# Add axis labels and title
plt.xlabel('Length of Time Series')
plt.ylabel('Error (%)')
plt.title('The Accuracy of DC Sparse DTW')

# Add legend
plt.legend()
file_name = 'DC_accuracy.png'
plt.savefig(file_name)




