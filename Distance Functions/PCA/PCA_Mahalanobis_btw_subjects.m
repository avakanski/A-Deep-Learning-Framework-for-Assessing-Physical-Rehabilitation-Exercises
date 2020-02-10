% PCA_Mahalanobis_btw_subjects: uses PCA to reduce the dimensionality of raw data, and 
% afterward calculates Mahalanobis distance on the low-dimensional data for the between-subject case

clear; clc; close all;

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

% Correct movement 
Train_Data_reduced = cell(1,T2);
for r =1:T2
    Train_Data_reduced{r}=Train_Data_rdim(:,(r-1)*T1+1:r*T1)';
end

% Incorrect movement
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

%% Scale and plot data

% Scale data
rms_train = rms_train/T1/T2;
rms_test = rms_test/T1/T2;

% Scale data in the [1,20] range
MAX = max(max(rms_train),max(rms_test));
MIN = min(min(rms_train),min(rms_test));
for s = 1: T2
    rms_train(s) = 19*(rms_train(s)-MIN)/(MAX-MIN)+1;
    rms_test(s) = 19*(rms_test(s)-MIN)/(MAX-MIN)+1;
end

% Plot
h = figure; 
plot(rms_train,'go','LineWidth',2); hold on, ...
plot(rms_test,'rs','LineWidth',2);
xlabel('Sequence Number', 'fontsize',18);
ylabel('Mahalanobis Distance', 'fontsize',18);
% Title: Mahalanobis Distance Between Subjects - PCA Dimensionality Reduction
title('Mahalanobis Distance BS PCA','fontsize',18);
legend({'Correct Sequences','Incorrect Sequences'}, 'fontsize',16,'location','NW')
set(gca,'box','off','fontweight','bold','LineWidth',2);
set(gcf,'Units','inches','position',[0 0 5.5 4.5]);
% print(h,'../../Results/Mahalanobis_Distance_BS_PCA','-dpng','-r300');

%% Separation degree betwween the correct and incorrect sequences
SD = zeros(1,1);
for i=1:T2
    for j=1:T2
        SD=SD+(rms_test(i)-rms_train(j))/(abs(rms_test(i))+abs(rms_train(j)));
    end
end
SD=SD/T2/T2;

