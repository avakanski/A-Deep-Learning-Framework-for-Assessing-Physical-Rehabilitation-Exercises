% Var_DTW_between_subjects: uses maximum variance to reduce the dimensionality of raw data, and 
% afterward calculates DTW distance on the low-dimensional data for the between-subject case

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

%% calculate the DTW between subjects
% correct sequences
dtw_train = zeros(1,T2);
for i=1:T2
    for j=1:T2
        dtw_train(i) =dtw_train(i)+dtw(Train_Var{i},Train_Var{j});
    end
end

% incorrect sequences
dtw_test = zeros(1,T2);  
for i=1:T2
    for j=1:T2
        dtw_test(i) =dtw_test(i)+dtw(Test_Var{i},Train_Var{j});
    end
end

%% scale data and plot it
dtw_train=dtw_train/T1/T2;
dtw_test=dtw_test/T1/T2;

% plot
h=figure;
plot(dtw_train,'bo');
hold on;
plot(dtw_test,'r*');
title('M1');
legend('Correct sequences','Incorrect sequences')
