% DTW_within_subjects: calculates DTW on the raw data for the within-subject case

clear; clc; close;

% Load the data for DeepSquat
load('../../Data for Distance Functions/M1-DeepSquat-Correct.mat');
load('../../Data for Distance Functions/M1-DeepSquat-Incorrect.mat');

% timesteps
T1 = size(Train_Data{1},1);
% repetitions
T2=length(Train_Data);

% repetions for each subject
rt=[9 9 9 9 9 10 8 9 8 10];

%% calculate the DTW
% correct sequences
dtw_train = zeros(1,T2);
M=0;
for s=1:10
    for i=1:rt(1,s)
        for j=1:rt(1,s)
            dtw_train(M+i)=dtw_train(M+i)+dtw(Train_Data{M+i},...
                               Train_Data{M+j});
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
            dtw_test(N+i) = dtw_test(N+i)+dtw(Test_Data{N+i},...
                               Train_Data{N+j});
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
