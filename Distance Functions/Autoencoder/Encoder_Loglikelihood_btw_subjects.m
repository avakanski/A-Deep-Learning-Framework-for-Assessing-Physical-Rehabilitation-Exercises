% Encoder_Loglikelihood_btw_subjects: uses autoencoder neural network to reduce the dimensionality of raw data, 
% and afterward calculates loglikelihood on the low-dimension data for the between-subject case

clear; clc; close;

%% Load the data

% Correct repetitions
Data_NN = csvread('../../Data/Autoencoder_Output_Correct.csv');

% Incorrect repetitions
Data_NN_inc = csvread('../../Data/Autoencoder_Output_Incorrect.csv');

% Data dimensions
nDim = 4;
% Number of timesteps
T1 = size(Data_NN,2)/nDim;
% Number of repetitions
T2 = size(Data_NN,1);

% Tranform the data into cells
% Correct repetitions
Train_Data_Reduced = cell(1,90);
for i=1:T2
    temp = [];
    for j=1:nDim
        temp = [temp; Data_NN(i,j:nDim:nDim*T1)];
    end
    Train_Data_Reduced{1,i} = temp';
end

% Incorrect repetitions
Test_Data_Reduced = cell(1,90);
for i=1:T2
    temp = [];
    for j=1:nDim
        temp = [temp; Data_NN_inc(i,j:nDim:nDim*T1)];
    end
    Test_Data_Reduced{1,i} = temp';
end

addpath('../../Utility Functions')

%% GMM encoding

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
% Title: GMM Loglikelihood Between Subjects - Autoencoder Dimensionality Reduction
title('GMM Loglikelihood BS ENC','fontsize',18);
legend({'Correct Sequences','Incorrect Sequences'}, 'fontsize',16,'location','NW')
set(gca,'box','off','fontweight','bold','LineWidth',2);
set(gcf,'Units','inches','position',[0 0 5.5 4.5]);
% print(h,'../../Results/GMM_Loglikelihood_BS_ENC','-dpng','-r300');

%% Separation degree betwween the correct and incorrect sequences

SD = zeros(1,1);
for i=1:T2
    for j=1:T2
        SD=SD+(loglikelihood_test(i)-loglikelihood_train(j))/(abs(loglikelihood_test(i))+abs(loglikelihood_train(j)));
    end
end
SD=SD/T2/T2;