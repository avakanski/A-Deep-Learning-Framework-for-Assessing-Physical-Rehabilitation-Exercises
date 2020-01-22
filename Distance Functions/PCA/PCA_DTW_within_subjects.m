% PCA_DTW_within_subjects: uses PCA to reduce the dimensionality of raw data, and 
% afterward calculates DTW distance on the low-dimensional data for the within-subject case

clear all; close all;

% Load the data for DeepSquat
load('../../Data for Distance Functions/M1-DeepSquat-Correct.mat');
load('../../Data for Distance Functions/M1-DeepSquat-Incorrect.mat');

% timesteps
T1 = size(Train_Data{1},1);
% repetitions
T2=length(Train_Data);

% repetions for each subject
rt=[9 9 9 9 9 10 8 9 8 10];

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

%% calculate the DTW
% correct sequences
dtw_train = zeros(1,T2);
M=0;
for s=1:10
    for i=1:rt(1,s)
        for j=1:rt(1,s)
            dtw_train(M+i)=dtw_train(M+i)+dtw(Train_Data_reduced{M+i},...
                               Train_Data_reduced{M+j});
        end
        % scale data
        if rt(1,s)>0
           dtw_train(M+i)=dtw_train(M+i)/rt(1,s);
        else
           dtw_train(M+i)=dtw_train(M+i);
        end
    end
    M=M+rt(1,s);
end

% incorrect sequences
dtw_test = zeros(1,T2);
N=0;
for s=1:10
    for i=1:rt(1,s)
        for j=1:rt(1,s)
            dtw_test(N+i) = dtw_test(N+i)+dtw(Test_Data_reduced{N+i},...
                               Train_Data_reduced{N+j});
        end
    % scale data
        if rt(1,s)>0
           dtw_test(N+i)=dtw_test(N+i)/rt(1,s);
        else
           dtw_test(N+i)=dtw_test(N+i);
        end
    end
    N=N+rt(1,s);
end

%% scale data and plot it
dtw_train=dtw_train/T1;
dtw_test=dtw_test/T1;

% plot
h=figure;
plot(dtw_train,'bo');
hold on;
plot(dtw_test,'r*');
title('M1');
legend('Correct sequences','Incorrect sequences')

