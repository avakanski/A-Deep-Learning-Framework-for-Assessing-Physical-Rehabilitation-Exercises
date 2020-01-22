% En_Euc_within_subjects: uses autoencoder neural network to reduce the dimensionality of raw data, 
% and afterward calculates Euclidean distance on the low-dimension data for the within-subject case

clear; clc; close all;

% Load the data (Part 1)
load('../../Data for Distance Functions/M1_Reduced-DeepSquat.mat');

% timesteps
T1 = size(Train_Data_Reduced{1},1);
% repetitions
T2=length(Train_Data_Reduced);

%repetitions for each subject
rt=[9 9 9 9 9 10 8 9 8 10];

%% Calculate RMS within subjects

% correct sequences
M=0;
rms_train=zeros(1,T2);
for n=1:10
   for i=1:rt(1,n)
       for j=1:rt(1,n)
           for t = 1:T1
               rms_train(M+i) = rms_train(M+i)+norm...
               (Train_Data_Reduced{M+i}(t,:)-Train_Data_Reduced{M+j}(t,:));
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
                (Test_Data_Reduced{N+i}(t,:)-Train_Data_Reduced{N+j}(t,:));
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
