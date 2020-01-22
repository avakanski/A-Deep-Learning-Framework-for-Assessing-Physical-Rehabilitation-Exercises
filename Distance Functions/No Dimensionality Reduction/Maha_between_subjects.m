% Maha_between_subjects: calculates Mahalanobis distance on the raw data for the between-subject case

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

%% M distance based on variances of each dimension

% Find the variances for each correct sequence
dim_var=zeros(T2,D);
for i = 1:T2
    dim_var(i,:) = var(Train_Data{i},1);
end


% Calculate the mean variance for all correct sequences
% to be used for calculating Mahalanobis distance

 Mhlnbs_dist = mean(dim_var);
 Mhlnbs_cov=diag(Mhlnbs_dist);
 

%% M distance for the correct sequences
rms_train = zeros(1,T2);
for i=1:T2 
    for j=1:T2
        for t = 1:T1
            rms_train(i) = rms_train(i)+sqrt((Train_Data{i}(t,:)...
            -Train_Data{j}(t,:))*Mhlnbs_cov*...
           (Train_Data{i}(t,:)-Train_Data{j}(t,:))');
        end
    end
end


%% Find the variances for each incorrect sequence 
dim_var_test=zeros(T2,D);
for i = 1:T2
    dim_var_test(i,:) = var(Test_Data{i},1);
end

%% M distance for the incorrect sequences
rms_test = zeros(1,T2);
for i=1:T2
    Mhlnbs_cov_test = diag(dim_var_test(i,:));
    for j=1:T2
        for t = 1:T1
            rms_test(i) = rms_test(i)+sqrt((Test_Data{i}(t,:)...
            -Train_Data{j}(t,:))*Mhlnbs_cov_test*...
            (Test_Data{i}(t,:)-Train_Data{j}(t,:))');
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
