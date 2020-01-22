% En_ DTW_within_subjects: uses autoencoder neural network to reduce the dimensionality of raw data, 
% and afterward calculates DTW distance on the low-dimensional data for the within-subject case

clear; close; clc;

% Load the data for DeepSquat
load('../../Data for Distance Functions/M1_Reduced-DeepSquat.mat');

% timesteps
T1 = size(Train_Data_Reduced{1},1);
% repetitions
T2=length(Train_Data_Reduced);

%% calculate the DTW
% correct sequences
dtw_train = zeros(1,T2);
for i=1:T2
    for j=1:T2
        dtw_train(i) =dtw_train(i)+dtw(Train_Data_Reduced{i},...
                      Train_Data_Reduced{j});
    end
end

% incorrect sequences
dtw_test = zeros(1,T2);  
for i=1:T2
    for j=1:T2
        dtw_test(i) =dtw_test(i)+dtw(Test_Data_Reduced{i},...
                      Train_Data_Reduced{j});
    end
end


%% scale data and plot it
% scale data
dtw_train=dtw_train/T1/T2;
dtw_test=dtw_test/T1/T2;

% plot
h=figure;
plot(dtw_train,'bo');
hold on;
plot(dtw_test,'r*');
title('M1');
legend('Correct sequences','Incorrect sequences')
