"""
This function is called to load data for the first exercise Deep Squat when RNN_GMM_Between_M1 is running.

"""

import csv
import numpy as np

# time steps
timesteps = 240
# repetition number
nr = 90

def load_data():
    f = open('../../Data for Neural Networks/M1_DeepSquat/Train_X1.csv')
    csv_f = csv.reader(f)
    Train_X = list(csv_f)

    # Convert the input sequences into numpy arrays
    train_input1 = np.asarray(Train_X, dtype = float)
    n_dim = 117
    train_input = np.zeros((nr,timesteps,n_dim))
    for i in range(len(train_input1)//n_dim):
          train_input[i,:,:] = np.transpose(train_input1[n_dim*i:n_dim*(i+1),:])
    
    f = open('../../Data for Neural Networks/M1_DeepSquat/Train_Y1.csv')
    csv_f = csv.reader(f)
    Train_Y = list(csv_f)
    
    # Convert the input labels into numpy arrays
    train_label = np.asarray(Train_Y, dtype = float)
    
    f = open('../../Data for Neural Networks/M1_DeepSquat/Test_X1.csv')
    csv_f = csv.reader(f)
    Test_X = list(csv_f)

    # Convert the input sequences into numpy arrays
    test_input1 = np.asarray(Test_X)
    n_dim = 117
    test_input = np.zeros((nr,timesteps,n_dim))
    for i in range(len(test_input1)//n_dim):
          test_input[i,:,:] = np.transpose(test_input1[n_dim*i:n_dim*(i+1),:])
            
    f = open('../../Data for Neural Networks/M1_DeepSquat/Test_Y1.csv')
    csv_f = csv.reader(f)
    Test_Y = list(csv_f)
    
    # Convert the input labels into numpy arrays
    test_label = np.asarray(Test_Y, dtype = float)
    
    return train_input, train_label, test_input, test_label