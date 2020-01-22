% Var_loglikelihood_between_subjects: uses maximum variance to reduce the dimensionality of raw data, and 
% afterward calculates loglikelihood on the low-dimensional data for the between-subject case

clear; clc; close;

% Load the data for DeepSquat
load('../../Data for Distance Functions/M1-DeepSquat-Correct.mat');
load('../../Data for Distance Functions/M1-DeepSquat-Incorrect.mat');

% timesteps
T1 = size(Train_Data{1},1);
% dimension
D = size(Train_Data{1},2);
% repetitions
T2=length(Train_Data);

%% Extract 3 dimensions with greatest variations
n_dim = 3;

% Find the variances for each correct sequence
dim_var=zeros(T2,D);
for i = 1:T2
    dim_var(i,:) = var(Train_Data{i},1);
end

% mean variance 
mean_var = mean(dim_var);

% Sort in descending order
bb = sort(mean_var,'descend');

% Extract the indices of the most varying dimensions
ind =zeros(1,n_dim);
i=1;
while i<=n_dim
    cc = find(mean_var == bb(i));
    % in case that two or above have the same variance
    for h=1:length(cc) 
         ind(1,i) = cc(h); 
    i = i+1;
    end
end
% cancell the redundant data that 
% when i=n_dim and length of cc is two or above
if length(ind)>n_dim
    ind = ind(1,1:n_dim);
end

% reduce correct sequences
Train_Var=cell(1,T2);
for i = 1:T2
    Train_Var{i} = Train_Data{i}(:,ind);
end

% reduce incorrect sequences
Test_Var=cell(1,T2);
for i = 1:T2
    Test_Var{i} = Test_Data{i}(:,ind);
end

%% get the distribution

% Create a row for the time indices of correct data
Data_Train=repmat(1:T1,1,T2);
Data_Train_position=[];
% Concatenate the data
for i=1:T2
    Data_Train_position=[Data_Train_position,Train_Var{i}'];
end
Data_Train=[Data_Train;Data_Train_position];


% Create a row for the time indices of incorrect data
Data_Test = repmat(1:T1,1,T2);
Data_Test_position=[];
% Concatenate the data
for i=1:T2
    Data_Test_position=[Data_Test_position,Test_Var{i}'];
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
