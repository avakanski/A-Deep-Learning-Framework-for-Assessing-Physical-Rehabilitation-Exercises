% PCA_Euc_between_subjects: uses PCA to reduce the dimensionality of raw data, and 
% afterward calculates Euclidean distance on the low-dimensional data for the between-subject case

clear; clc; close all;

% Load the data (Part 1)
load('../../Data for Distance Functions/M1-DeepSquat-Correct.mat');
load('../../Data for Distance Functions/M1-DeepSquat-Incorrect.mat');

% timesteps
T1 = size(Train_Data{1},1);
% repetitions
T2=length(Train_Data);


%% Reshape the data Part 2)

% Concatenate the data of correct movements
Train_DataR=[];
for i=1:T2
    Train_DataR=[Train_DataR,Train_Data{i}'];
end

% Concatenate the data of incorrect movements
Test_DataR=[];
for i=1:T2
    Test_DataR=[Test_DataR,Test_Data{i}'];
end


%% Perform PCA (Principal Component Analysis) (Part 3)

% Number of principal components.
nbPC = 3;

% Extract the eigenvectors E and eigenvalues v of the covariance matrix
[E,v] = eig(cov(Train_DataR'));
E = fliplr(E);
% Compute the transformation matrix by keeping the first nbPC eigenvectors
A = E(:,1:nbPC);
% Project the training data in the latent space
Train_Data_rdim = A' * Train_DataR;
% Project the test data in the latent space
Test_Data_rdim = A' * Test_DataR;

%% recover data structure
% correct sequences 
Train_Data_reduced = cell(1,T2);
for r =1:T2
    Train_Data_reduced{r}=Train_Data_rdim(:,(r-1)*T1+1:r*T1)';
end

% incorrect sequences
Test_Data_reduced = cell(1,T2);
for r =1:T2
    Test_Data_reduced{r}=Test_Data_rdim(:,(r-1)*T1+1:r*T1)';
end

%% RMS between subjects
% correct sequences
rms_train=zeros(1,T2);
for i=1:T2
    for j=1:T2
        for t = 1:T1
            rms_train(i) = rms_train(i)+norm...
            (Train_Data_reduced{i}(t,:)-Train_Data_reduced{j}(t,:));
        end
    end
end
   
% incorrect sequences
rms_test=zeros(1,T2);  
for i=1:T2
    for j=1:T2
        for t = 1:T1
            rms_test(i) = rms_test(i)+norm...
            (Train_Data_reduced{i}(t,:)-Test_Data_reduced{j}(t,:));
        end
    end
end
  

%% scale data and plot it
rms_train = rms_train/T1/T2;
rms_test = rms_test/T1/T2;

% plot
h=figure;
plot(rms_train,'bo');
hold on;
plot(rms_test,'r*');
title('M1');
legend('Correct sequences','Incorrect sequences')
