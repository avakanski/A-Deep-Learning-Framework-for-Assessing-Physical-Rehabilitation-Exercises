% PCA_loglikelihood_btw_subjects: uses PCA to reduce the dimensionality of raw data, and 
% afterward calculates loglikelihood on the low-dimensional data for the between-subject case

clear; clc; close;

%% Load the data 

% Correct repetitions
Data_NN = csvread('../../Data/Data_Correct.csv');

% Incorrect repetitions
Data_NN_inc = csvread('../../Data/Data_Incorrect.csv');

% Data dimensions
n_dim = 117;
% Number of timesteps
T1 = size(Data_NN,2);
% Number of repetitions
T2 = size(Data_NN,1)/n_dim;

% Tranform the data into cells
% Correct repetitions
Train_Data = cell(1,90);
for i=1:T2
    Train_Data{1,i} = Data_NN((i-1)*n_dim+1:i*n_dim,:)';
end

% Incorrect repetitions
Test_Data = cell(1,90);
for i=1:T2
    Test_Data{1,i} = Data_NN_inc((i-1)*n_dim+1:i*n_dim,:)';
end

addpath('../../Utility Functions')

%% Reshape the data

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


%% Perform PCA (Principal Component Analysis)

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

%% Recover data structure

% Correct movements 
Train_Data_reduced = cell(1,T2);
for r =1:T2
    Train_Data_reduced{r}=Train_Data_rdim(:,(r-1)*T1+1:r*T1)';
end

% Incorrect movements
Test_Data_reduced = cell(1,T2);
for r =1:T2
    Test_Data_reduced{r}=Test_Data_rdim(:,(r-1)*T1+1:r*T1)';
end

%% GMM encoding

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
[Priors, Mu, Sigma] = EM_init_regularTiming(Data_Train,nbStates);
[Priors, Mu, Sigma] = EM_boundingCov(Data_Train,Priors, Mu, Sigma);

disp(['GMM modeling has been completed!',char(10)]);

%% Data loglikelihood

% Correct sequences
loglikelihood_train=zeros(1,T2);
for j=1:T2
    loglikelihood_train(j) = -loglik(Data_Train(:,(j-1)*T1...
    +1:j*T1), nbStates, Priors, Mu, Sigma);
end

% Incorrect sequences
loglikelihood_test=zeros(1,T2);
for j=1:T2
    loglikelihood_test(j) = -loglik(Data_Test(:,(j-1)*T1...
    +1:j*T1), nbStates, Priors, Mu, Sigma);
end

%% Scale and plot data

% don't need to divide loglikelihood_train by T1*T2
% because we have done this when calculate log-likelihood distance

% Scale data in the [1,20] range
MAX = max(max(loglikelihood_train),max(loglikelihood_test));
MIN = min(min(loglikelihood_train),min(loglikelihood_test));
for s = 1: T2
    loglikelihood_train(s) = 19*(loglikelihood_train(s)-MIN)/(MAX-MIN)+1;
    loglikelihood_test(s) = 19*(loglikelihood_test(s)-MIN)/(MAX-MIN)+1;
end

% Plot
h = figure; 
plot(loglikelihood_train,'go','LineWidth',2); hold on, ...
plot(loglikelihood_test,'rs','LineWidth',2);
xlabel('Sequence Number', 'fontsize',18);
ylabel('GMM Loglikelihood', 'fontsize',18);
% Title: GMM Loglikelihood Between Subjects - PCA Dimensionality Reduction
title('GMM Loglikelihood BS PCA','fontsize',18);
legend({'Correct Sequences','Incorrect Sequences'}, 'fontsize',16,'location','NW')
set(gca,'box','off','fontweight','bold','LineWidth',2);
set(gcf,'Units','inches','position',[0 0 5.5 4.5]);
% print(h,'../../Results/GMM_Loglikelihood_BS_PCA','-dpng','-r300');

%% Separation degree betwween the correct and incorrect sequences
SD = zeros(1,1);
for i=1:T2
    for j=1:T2
        SD=SD+(loglikelihood_test(i)-loglikelihood_train(j))/(abs(loglikelihood_test(i))+abs(loglikelihood_train(j)));
    end
end
SD=SD/T2/T2;
