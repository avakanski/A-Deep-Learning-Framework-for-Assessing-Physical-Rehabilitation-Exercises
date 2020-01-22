% Var_Maha_between_subjects: uses maximum variance to reduce the dimensionality of raw data, and 
% afterward calculates Mahalanobis distance on the low-dimensional data for the between-subject case

clear; clc; close;

% Load the data
load('../../Data for Distance Functions/M1-DeepSquat-Correct.mat');
load('../../Data for Distance Functions/M1-DeepSquat-Incorrect.mat');

% timesteps
T1 = size(Train_Data{1},1);
% dimension
D = size(Train_Data{1},2);
% repetitions
T2=length(Train_Data);

%% Extract 3 dimensions with greatest variations(reduce data)
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

%% Mahalanobis distance based on variances of each dimension

% Find the variances for each correct sequence(reduced)
dim_var=zeros(T2,n_dim);
for i = 1:T2
    dim_var(i,:) = var(Train_Var{i},1);
end

% Calculate the mean variance for all correct sequences
% to be used for calculating Mahalanobis distance
Mhlnbs_dist = mean(dim_var);
Mhlnbs_cov=diag(Mhlnbs_dist);

% M distance for the correct sequence
rms_train = zeros(1,T2);
for i=1:T2 
    for j=1:T2
        for t = 1:T1
            rms_train(i) = rms_train(i)+sqrt((Train_Var{i}(t,:)...
            -Train_Var{j}(t,:))*Mhlnbs_cov*...
           (Train_Var{i}(t,:)-Train_Var{j}(t,:))');
        end
    end
end


%% Find the variances for each incorrect sequence 
dim_var_test=zeros(T2,n_dim);
for i = 1:T2
    dim_var_test(i,:) = var(Test_Var{i},1);
end

% M distance for the incorrect sequence
rms_test = zeros(1,T2);
for i=1:T2
    Mhlnbs_cov_test = diag(dim_var_test(i,:));
    for j=1:T2
        for t = 1:T1
            rms_test(i) = rms_test(i)+sqrt((Test_Var{i}(t,:)...
            -Train_Var{j}(t,:))*Mhlnbs_cov_test*...
            (Test_Var{i}(t,:)-Train_Var{j}(t,:))');
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
