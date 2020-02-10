% PCA_loglikelihood_wthn_subjects: uses PCA to reduce the dimensionality of raw data, and 
% afterward calculates loglikelihood on the low-dimensional data for the within-subject case

clear; close; clc;

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

% repetion for subject
rt=[9 9 9 9 9 10 8 9 8 10];

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
Data_Train=cell(1,10);
Data_Train_position=cell(1,10);
M=0;
for s=1:10   
    Data_Train{s}=repmat(1:T1,1,rt(1,s));
    % Concatenate the data
    for i=1:rt(1,s)
        Data_Train_position{s}=[Data_Train_position{s},Train_Data_reduced{M+i}'];
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
        Data_Test_position{s}=[Data_Test_position{s},Test_Data_reduced{N+i}'];
    end
    Data_Test{s} = [Data_Test{s};Data_Test_position{s}];
    N=N+rt(1,s);
end

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
       [Priors{1,s}, Mu{1,s}, Sigma{1,s}] = EM_init_regularTiming...
       (Data_Train{s},nbStates);
       [Priors{s}, Mu{s}, Sigma{s}] = EM_boundingCov(Data_Train{s},...
                               Priors{s}, Mu{s}, Sigma{s});
    end
end

disp(['GMM modeling has been completed!',char(10)]);

%% Data loglikelihood

% Correct sequences
loglikelihood_train=zeros(1,T2);
M=0;
for s=1:10
    for j=1:rt(1,s)
        loglikelihood_train(M+j) = -loglik(Data_Train{s}(:,...
        (j-1)*T1+1:j*T1),...
        nbStates, Priors{1,s}, Mu{1,s}, Sigma{1,s});
    end
    M=M+rt(1,s);
end

% Incorrect sequences
loglikelihood_test=zeros(1,T2);
N=0;
for s=1:10
    for j=1:rt(1,s)
        loglikelihood_test(N+j) = -loglik(Data_Test{s}(:,...
        (j-1)*T1+1:j*T1),...
         nbStates, Priors{1,s}, Mu{1,s}, Sigma{1,s});
    end
    N=N+rt(1,s);
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
% Title: GMM Loglikelihood Within Subjects - PCA Dimensionality Reduction
title('GMM Loglikelihood WS PCA','fontsize',18);
legend({'Correct Sequences','Incorrect Sequences'}, 'fontsize',16,'location','NW')
set(gca,'box','off','fontweight','bold','LineWidth',2);
set(gcf,'Units','inches','position',[0 0 5.5 4.5]);
% print(h,'../../Results/GMM_Loglikelihood_WS_PCA','-dpng','-r300');

%% Separation degree betwween the correct and incorrect sequences
SD = zeros(1,1);
for i=1:T2
    for j=1:T2
        SD=SD+(loglikelihood_test(i)-loglikelihood_train(j))/(abs(loglikelihood_test(i))+abs(loglikelihood_train(j)));
    end
end
SD=SD/T2/T2;

