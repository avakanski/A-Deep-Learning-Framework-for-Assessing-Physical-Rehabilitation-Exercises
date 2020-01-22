"""
This function is called to load the augmented data when RNN_GMM_Between_M1_Aug is running.

"""

import csv
import numpy as np

timesteps = 240

# repetition number
nr = 90

def load_data(a):
  
    f = open('../../Data for Neural Networks/M1_Aug_DeepSquat/X' + str(a)+'_movement1.csv')
    csv_f = csv.reader(f)
    Train_X = list(csv_f)

    # Convert the input sequences into numpy arrays
    train_input1 = np.asarray(Train_X, dtype = float)
    n_dim = 117
    train_input = np.zeros((nr,timesteps,n_dim))
    for i in range(len(train_input1)//n_dim):
          train_input[i,:,:] = np.transpose(train_input1[n_dim*i:n_dim*(i+1),:])
            
    
    ############################################################################ 
    
    f = open('../../Data for Neural Networks/M1_Aug_DeepSquat/Y'+str(a)+'_movement1.csv')
    csv_f = csv.reader(f)
    Train_Y = list(csv_f)
    
    # Convert the input labels into numpy arrays
    train_label = np.asarray(Train_Y, dtype = float)
    
    return train_input, train_label