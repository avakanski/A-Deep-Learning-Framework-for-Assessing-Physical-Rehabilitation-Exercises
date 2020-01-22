% Var_loglikelihood_within_subjects: uses maximum variance to reduce the dimensionality of raw data, and 
% afterward calculates loglikelihood on the low-dimensional data for the within-subject case

clear; close; clc;

% Load the data 
load('../../Data for Distance Functions/M1-DeepSquat-Correct.mat');
load('../../Data for Distance Functions/M1-DeepSquat-Incorrect.mat');

% timesteps
T1 = size(Train_Data{1},1);
% dimension
D = size(Train_Data{1},2);
% repetitions
T2=length(Train_Data);
% repetion for subject
rt=[9 9 9 9 9 10 8 9 8 10];

%% Extract 3 dimensions with greatest variations (Part 2)
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

%% Create a row for the time indices of correct data
Data_Train=cell(1,10);
Data_Train_position=cell(1,10);
M=0;
for s=1:10   
    Data_Train{s}=repmat(1:T1,1,rt(1,s));
    % Concatenate the data
    for i=1:rt(1,s)
        Data_Train_position{s}=[Data_Train_position{s},Train_Var{M+i}'];
    end
    Data_Train{s}=[Data_Train{s};Data_Train_position{s}];
    M=M+rt(1,s);
end

% Create a row for the time indices of incorrect data
Data_Test=cell(1,10);
Data_Test_position=cell(1,10);
N=0;
for s=1:10
    Data_Test{s} = repmat(1:T1,1,rt(1,s));
    % Concatenate the data
    for i=1:rt(1,s)
        Data_Test_position{s}=[Data_Test_position{s},Test_Var{N+i}'];
    end
    Data_Test{s} = [Data_Test{s};Data_Test_position{s}];
    N=N+rt(1,s);
end


%% Train GMM model
% Define the number of states for GMM
nbStates = 6;
nbVar=zeros(1,10);
Priors=cell(1,10);
Mu=cell(1,10);
Sigma=cell(1,10);

for s=1:10
    nbVar(1,s) = size(Data_Train{s},1);
    if nbVar(1,s)>1
       % Training by EM algorithm, initialized by k-means clustering.
       % [Priors, Mu, Sigma] = EM_init_kmeans(Data, nbStates);
       [Priors{1,s}, Mu{1,s}, Sigma{1,s}] = EM_init_regularTiming...
       (Data_Train{s},nbStates);
       % [Priors, Mu, Sigma] = EM(Data, Priors, Mu, Sigma);
       [Priors{s}, Mu{s}, Sigma{s}] = EM_boundingCov(Data_Train{s},...
                               Priors{s}, Mu{s}, Sigma{s});
    end
end

disp(['GMM modeling has been completed!',char(10)]);



%% Calculate loglikelihood (Part 4)
% correct sequences
loglikelihood_train=zeros(1,T2);
M=0;
for s=1:10
    for j=1:rt(1,s)
        loglikelihood_train(M+j) = loglik(Data_Train{s}(:,...
        (j-1)*T1+1:j*T1),...
        nbStates, Priors{1,s}, Mu{1,s}, Sigma{1,s});
    end
    M=M+rt(1,s);
end

% incorrect sequences
loglikelihood_test=zeros(1,T2);
N=0;
for s=1:10
    for j=1:rt(1,s)
        loglikelihood_test(N+j) = loglik(Data_Test{s}(:,...
        (j-1)*T1+1:j*T1),...
         nbStates, Priors{1,s}, Mu{1,s}, Sigma{1,s});
    end
    N=N+rt(1,s);
end

%% scale data and plot it
loglikelihood_train = -loglikelihood_train;
loglikelihood_test = -loglikelihood_test;

% plot
h=figure;
plot(loglikelihood_train,'bo');
hold on;
plot(loglikelihood_test,'r*');
title('M1');
legend('Correct sequences','Incorrect sequences')
