% PCA_loglikelihood_between_subjects: uses PCA to reduce the dimensionality of raw data, and 
% afterward calculates loglikelihood on the low-dimensional data for the between-subject case

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
% correct movement 
Train_Data_reduced = cell(1,T2);
for r =1:T2
    Train_Data_reduced{r}=Train_Data_rdim(:,(r-1)*T1+1:r*T1)';
end

% incorrect movement
Test_Data_reduced = cell(1,T2);
for r =1:T2
    Test_Data_reduced{r}=Test_Data_rdim(:,(r-1)*T1+1:r*T1)';
end

%% get the distribution
% Create a row for the time indices of correct data
Data_Train=repmat(1:T1,1,T2);
Data_Train_position=[];
% Concatenate the data
for i=1:T2
    Data_Train_position=[Data_Train_position,Train_Data_reduced{i}'];
end
Data_Train=[Data_Train;Data_Train_position];


% Create a row for the time indices of incorrect data
Data_Test = repmat(1:T1,1,T2);
Data_Test_position=[];
% Concatenate the data
for i=1:T2
    Data_Test_position=[Data_Test_position,Test_Data_reduced{i}'];
end
Data_Test = [Data_Test;Data_Test_position]; 

% Define the number of states for GMM
nbStates = 6;

% Train GMM
nbVar = size(Data_Train,1);
% Training by EM algorithm, initialized by k-means clustering.
% [Priors, Mu, Sigma] = EM_init_kmeans(Data, nbStates);
[Priors, Mu, Sigma] = EM_init_regularTiming(Data_Train,nbStates);
% [Priors, Mu, Sigma] = EM(Data, Priors, Mu, Sigma);
[Priors, Mu, Sigma] = EM_boundingCov(Data_Train,Priors, Mu, Sigma);

disp(['GMM modeling has been completed!',char(10)]);

%% Data loglikelihood (Part 4)
% correct sequences
loglikelihood_train=zeros(1,T2);
for j=1:T2
    loglikelihood_train(j) = loglik(Data_Train(:,(j-1)*T1...
    +1:j*T1), nbStates, Priors, Mu, Sigma);
end

% incorrect sequences
loglikelihood_test=zeros(1,T2);
for j=1:T2
    loglikelihood_test(j) = loglik(Data_Test(:,(j-1)*T1...
    +1:j*T1), nbStates, Priors, Mu, Sigma);
end

%% scale data and plot it
% don't need to divide loglikelihood_train{m} by T1(m)*T2(m)
% because we have done this when calculate log-likelihood distance
loglikelihood_train = -loglikelihood_train;
loglikelihood_test = -loglikelihood_test;

% plot
h=figure;
plot(loglikelihood_train,'bo');
hold on;
plot(loglikelihood_test,'r*');
title('M1');
legend('Correct sequences','Incorrect sequences')
