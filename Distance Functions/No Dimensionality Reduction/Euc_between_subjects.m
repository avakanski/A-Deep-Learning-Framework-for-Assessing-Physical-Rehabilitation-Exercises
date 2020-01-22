% Euc_between_subjects: calculates Euclidean distance on the raw data for the between-subject case

clear;clc;

% load data
load('../../Data for Distance Functions/M1-DeepSquat-Correct.mat');
load('../../Data for Distance Functions/M1-DeepSquat-Incorrect.mat');

% timesteps
T1 = size(Train_Data{1},1);
% repetitions
T2=length(Train_Data);


%% calculate RMS

% correct sequences
rms_train = zeros(1,T2);
for i=1:T2
    for j=1:T2
    for t = 1:T1
        rms_train(i) = rms_train(i)+norm(Train_Data{i}(t,:)...
                            -Train_Data{j}(t,:));
    end
    end
end
  
% incorrect sequences  
rms_test=zeros(1,T2);
for i=1:T2
    for j=1:T2
    for t = 1:T1
        rms_test(i) = rms_test(i)+norm(Test_Data{i}(t,:)...
                           -Train_Data{j}(t,:));
    end
    end
end


%% scale and plot data

% scale data
rms_train = rms_train/T1/T2;
rms_test = rms_test/T1/T2;

% plot
h=figure;
plot(rms_train,'bo');
hold on;
plot(rms_test,'r*');
title('M1');
legend('Correct sequences','Incorrect sequences')
