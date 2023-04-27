import os
import pandas as pd
import scipy
from scipy import io

base_path = 'E:/graduation project/data/ECG/archive/Training_2'
def mat_to_csv(file):
    '''
    A tool to concert a mat file to a csv file
    :param file: File name
    :return: data frame of the converted csv file
    '''
    features_struct = scipy.io.loadmat(file)
    features = list(features_struct.values())[-1]
    dfdata = pd.DataFrame(data=features)
    return dfdata