% PCA_Euc_within_subjects: uses PCA to reduce the dimensionality of raw data, and 
% afterward calculates Euclidean distance on the low-dimensional data for the within-subject case

clear; clc; close all;

% Load the data (Part 1)
load('../../Data for Distance Functions/M1-DeepSquat-Correct.mat');
load('../../Data for Distance Functions/M1-DeepSquat-Incorrect.mat');

% timesteps
T1 = size(Train_Data{1},1);
% repetitions
T2=length(Train_Data);

% repetition for each subject
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

%% Calculate RMS within subjects
% correct sequences
M=0;
rms_train=zeros(1,T2);
for n=1:10
   for i=1:rt(1,n)
       for j=1:rt(1,n)
           for t = 1:T1
               rms_train(M+i) = rms_train(M+i)+norm...
               (Train_Data_reduced{M+i}(t,:)-Train_Data_reduced{M+j}(t,:));
           end
       end
        % scale data
        if rt(1,n)>0
           rms_train(M+i)=rms_train(M+i)/rt(1,n);
        else
           rms_train(M+i)=rms_train(M+i);
        end
    end
    M=M+rt(1,n);
end

% incorrect sequences
N=0;
rms_test=zeros(1,T2);
for n=1:10
    for i=1:rt(1,n)
        for j=1:rt(1,n)
            for t = 1:T1
                rms_test(N+i)= rms_test(N+i)+norm...
                (Test_Data_reduced{N+i}(t,:)-Train_Data_reduced{N+j}(t,:));
            end   
        end
        % scale data
        if rt(1,n)>0
           rms_test(N+i)=rms_test(N+i)/rt(1,n);
        else
           rms_test(N+i)=rms_test(N+i);
        end
    end
    N=N+rt(1,n);
end

%% scale data and draw pictures
% scale data
rms_train = rms_train/T1;
rms_test = rms_test/T1;

% plot
h=figure;
plot(rms_train,'bo');
hold on;
plot(rms_test,'r*');
title('M1');
legend('Correct sequences','Incorrect sequences')

