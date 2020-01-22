% PCA_Maha_between_subjects: uses PCA to reduce the dimensionality of raw data, and 
% afterward calculates Mahalanobis distance on the low-dimensional data for the between-subject case

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


%% Mahalanobis distance based on variances of each dimension

% Find the variances for each correct sequence
dim_var=zeros(T2,nbPC);
for i = 1:T2
    dim_var(i,:) = var(Train_Data_reduced{i},1);
end

% Calculate the mean variance for all correct sequences
% to be used for calculating Mahalanobis distance
 Mhlnbs_dist = mean(dim_var);
 Mhlnbs_cov=diag(Mhlnbs_dist);

% M distance for the correct sequences
rms_train = zeros(1,T2);
for i=1:T2 
    for j=1:T2
        for t = 1:T1
            rms_train(i) = rms_train(i)+sqrt((Train_Data_reduced{i}(t,:)...
            -Train_Data_reduced{j}(t,:))*Mhlnbs_cov*...
           (Train_Data_reduced{i}(t,:)-Train_Data_reduced{j}(t,:))');
        end
    end
end


%% Find the variances for each incorrect sequence 
dim_var_test=zeros(T2,nbPC);
for i = 1:T2
    dim_var_test(i,:) = var(Test_Data_reduced{i},1);
end

% M distance for the incorrect sequences
rms_test = zeros(1,T2);
for i=1:T2
    Mhlnbs_cov_test = diag(dim_var_test(i,:));
    for j=1:T2
        for t = 1:T1
            rms_test(i) = rms_test(i)+sqrt((Test_Data_reduced{i}(t,:)...
            -Train_Data_reduced{j}(t,:))*Mhlnbs_cov_test*...
            (Test_Data_reduced{i}(t,:)-Train_Data_reduced{j}(t,:))');
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

