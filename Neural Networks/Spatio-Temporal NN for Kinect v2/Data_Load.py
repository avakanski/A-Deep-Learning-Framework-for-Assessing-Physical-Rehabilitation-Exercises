"""
The file loads the data for full body skeletons

"""

import csv
import numpy as np

def load_data():
    f = open('Data_KIMORE_e5/Train_X.csv')
    csv_f = csv.reader(f)
    Train_X = list(csv_f)

    # Convert the input sequences into numpy arrays
    train_input1 = np.asarray(Train_X, dtype = float)
    n_dim = 88
    train_input = np.zeros((204,100,n_dim))
    for i in range(len(train_input1)//100):
          train_input[i,:,:] = train_input1[100*i:100*(i+1),:]
    
    f = open('Data_KIMORE_e5/Train_Y.csv')
    csv_f = csv.reader(f)
    Train_Y = list(csv_f)
    
    # Convert the input labels into numpy arrays
    train_label = np.transpose(np.asarray(Train_Y[0], dtype = float))
    
    return train_input, train_label