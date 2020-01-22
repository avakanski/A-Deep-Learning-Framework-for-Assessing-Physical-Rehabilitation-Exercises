% En_Loglikelihood_between_subjects: uses autoencoder neural network to reduce the dimensionality of raw data, 
% and afterward calculates loglikelihood on the low-dimension data for the between-subject case

clear; clc; close;

% Load the data for DeepSquat
load('../../Data for Distance Functions/M1_Reduced-DeepSquat.mat');

% timesteps
T1 = size(Train_Data_Reduced{1},1);
% repetitions
T2=length(Train_Data_Reduced);

%% get the distribution

% Create a row for the time indices of correct data
Data_Train=repmat(1:T1,1,T2);
Data_Train_position=[];
% Concatenate the data
for i=1:T2
    Data_Train_position=[Data_Train_position,Train_Data_Reduced{i}'];
end
Data_Train=[Data_Train;Data_Train_position];


% Create a row for the time indices of incorrect data
Data_Test = repmat(1:T1,1,T2);
Data_Test_position=[];
% Concatenate the data
for i=1:T2
    Data_Test_position=[Data_Test_position,Test_Data_Reduced{i}'];
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
